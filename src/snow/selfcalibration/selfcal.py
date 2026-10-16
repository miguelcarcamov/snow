from __future__ import annotations

import copy
import os
import shutil
import warnings
from abc import ABCMeta, abstractmethod
from pathlib import Path
from dataclasses import dataclass

from casatasks import (
    applycal, clearcal, delmod, flagdata, flagmanager, gaincal, split, statwt, uvsub
)

from ..imaging.imager import Imager
from ..utils.selfcal_utils import (
    columns_match, count_spectral_windows, get_spw_ids, is_column_in_ms
)
from .caltable import CalTable, build_gaintable_arguments, derive_spwmap

# Data selection shared by the solve, the apply and the imager (#27)
SELECTION_PARAMETERS = ("field", "spw")
# flagdata modes that flag outliers in the residual visibilities
FLAG_MODES = ("rflag", "tfcrop")


@dataclass(init=False, repr=True)
class Selfcal(metaclass=ABCMeta):

    def __init__(
        self,
        visfile: str = "",
        imager: Imager = None,
        field: str = "",
        spw: str = "",
        refant: str = "",
        spwmap: list = [],
        minblperant: int = 4,
        want_plot: bool = True,
        interp: str = 'linear',
        gaintype: str = 'T',
        uvrange: str = "",
        solint: list = [],
        varchange_imager: dict = None,
        varchange_selfcal: dict = None,
        output_caltables: str = None,
        previous_selfcal: Selfcal = None,
        input_caltable: str = None,
        minsnr: float = 3.0,
        applymode: str = "calflag",
        flag_mode: str = "rflag",
        combine: str = "",
        flag_dataset: bool = False,
        restore_psnr: bool = False,
        subtract_source: bool = False
    ):
        """
        General self-calibration class

        Parameters
        ----------
        visfile :
            Input visibility measurement set
        imager :
            Imager instance to run. It must write the model to MODEL_DATA; save_model is forced on
        field :
            Field selection used to solve, apply and image. An empty selection takes the imager's
        spw :
            Spectral window selection used to solve, apply and image. An empty selection takes the
            imager's
        refant :
            Reference antenna
        spwmap :
            Spectral window map
        minblperant :
            Minimum baseline per antenna
        want_plot :
            Whether to plot the calibration tables or not
        interp :
            Temporal interpolation for each gaintable
        gaintype :
            Type of gain solution ("G", "T", or "GSPLINE")
        uvrange :
            Select data within uvrange (default units meters)
        solint :
            Solution interval: e.g "inf", "60s", "int"
        varchange_imager :
            Dictionary of imager variables that change on each iteration
        varchange_selfcal :
            Dictionary of self-cal variables that change on each iteration
        output_caltables :
            Output path to save calibration tables
        previous_selfcal :
            Previous self-cal object if any
        input_caltable :
            Input calibration table if any
        minsnr :
            Reject solutions below this SNR
        applymode :
            Calibration mode - ””=”calflag”, ”calflagstrict”, ”trial”, ”flagonly”, ”flagonlystrict”, or ”calonly”.
            "calflag" flags the data whose solutions fail minsnr, and these flags add up over
            iterations; "calonly" keeps those data uncalibrated but unflagged
        flag_mode :
            Mode used to flag residual outliers: "rflag" or "tfcrop"
        combine :
            Data axes to combine for solving
        flag_dataset :
            Whether to flag Fourier residuals outliers
        restore_psnr :
            Restores the dataset if the peak signal-to-noise ratio decreases
        subtract_source :
            Subtract source model if needed
        """
        # Continuing from another self-calibration must start from its measurement set; any other
        # dataset silently discards its solutions and model (#25).
        if previous_selfcal is not None:
            if visfile == "":
                visfile = previous_selfcal.visfile
            elif Path(visfile).resolve() != Path(previous_selfcal.visfile).resolve():
                raise ValueError(
                    "visfile ({0}) is not the measurement set of previous_selfcal ({1}). Starting "
                    "from another dataset discards the previous self-calibration. Pass "
                    "visfile=previous_selfcal.visfile (or leave visfile empty), and create this "
                    "object after previous_selfcal.run().".format(
                        visfile, previous_selfcal.visfile
                    )
                )

        # Public variables
        self.visfile = visfile
        self.imager = imager
        self.field = field
        self.spw = spw
        self.refant = refant
        self.spwmap = spwmap
        self.minblperant = minblperant
        self.want_plot = want_plot
        self.interp = interp
        self.gaintype = gaintype
        self.uvrange = uvrange
        self.solint = solint
        self.varchange_imager = varchange_imager
        self.varchange_selfcal = varchange_selfcal
        self.output_caltables = output_caltables
        self.previous_selfcal = previous_selfcal
        self.input_caltable = input_caltable
        self.minsnr = minsnr
        self.applymode = applymode
        self.flag_mode = flag_mode
        self.combine = combine
        self.flag_dataset = flag_dataset
        self.restore_psnr = restore_psnr
        self.subtract_source = subtract_source

        # Protected variables
        self._caltables = []
        # Accepted tables currently applied to the data, in application order (#26)
        self._caltable_stack = []
        self._caltables_versions = []
        self._psnr_history = []
        # Fraction of the selected data flagged before self-calibration and after each iteration
        self._flagged_history = []
        self._calmode = ""
        self._loops = 0
        self._psnr_visfile_backup = self.visfile

        if self.imager is None:
            self._image_name = ""

        if output_caltables is None:
            self.output_caltables = self.imager.output

        if self.varchange_imager is not None:
            list_of_values = [value for key, value in self.varchange_imager.items()]
            it = iter(list_of_values)
            if not all(len(l) == len(self.solint) for l in it):
                raise ValueError(
                    "Error, length of solint and variable that changes through iterations must be the same"
                )

        if self.varchange_selfcal is not None:
            list_of_values = [value for key, value in self.varchange_selfcal.items()]
            it = iter(list_of_values)
            if not all(len(l) == len(self.solint) for l in it):
                raise ValueError(
                    "Error, length of solint and variable that changes through iterations must be the same"
                )

        self._sync_selection()

        if self.flag_dataset and self.flag_mode not in FLAG_MODES:
            raise ValueError(
                "flag_mode '{0}' cannot flag residual outliers; use one of {1}".format(
                    self.flag_mode, FLAG_MODES
                )
            )

        if self.subtract_source:
            if self.imager.getPhaseCenter() != "":
                raise ValueError(
                    "Error, phase center needs to be set if a source is going to be subtracted"
                )

    @property
    def imager(self):
        return self.__imager

    @imager.setter
    def imager(self, input_imager):
        if input_imager is not None:
            if not isinstance(input_imager, Imager):
                raise ValueError("The imager attribute must be an instance of Imager")
            else:
                self.__imager = input_imager
                self.__imager.inputvis = self.visfile
                # Without a model, gaincal solves against a 1 Jy point source at the phase centre
                if not input_imager.save_model:
                    warnings.warn(
                        "Self-calibration needs the model in MODEL_DATA: setting save_model=True "
                        "on the imager."
                    )
                    self.__imager.save_model = True
                self._image_name = self.__imager.output
        else:
            self.__imager = None

    @property
    def input_caltable(self):
        return self.__input_caltable

    @input_caltable.setter
    def input_caltable(self, caltable):
        if caltable is not None:
            if isinstance(caltable, str):
                if caltable != "":
                    if not os.path.exists(caltable):
                        raise FileNotFoundError(
                            "The caltable " + self.input_caltable + " needs to be created"
                        )
                    else:
                        self.__input_caltable = caltable
                else:
                    self.__input_caltable = caltable
            else:
                raise ValueError("The input_caltable should be a string")
        else:
            self.__input_caltable = ""

    def _sync_selection(self) -> None:
        """
        Protected method that makes solving, applying and imaging use one field and spectral window
        selection (#27). Otherwise part of the data is calibrated without matching solutions, or
        solutions are found from data that are never imaged. An empty selection takes the other.
        """
        for name in SELECTION_PARAMETERS:
            own, imager_value = getattr(self, name), getattr(self.imager, name)
            if own and imager_value and own != imager_value:
                raise ValueError(
                    "{0}='{1}' differs from the imager {0}='{2}'. Solving, applying and imaging "
                    "must use the same selection: set it in one place.".format(
                        name, own, imager_value
                    )
                )
            setattr(self, name, own or imager_value)
            setattr(self.imager, name, own or imager_value)

    def _warn_if_corrected_differs_from_data(self) -> None:
        """
        Protected method that warns when the input measurement set already holds a calibration in
        CORRECTED_DATA. The imager images CORRECTED_DATA, but gaincal solves against DATA without
        that calibration, so the input must be split to put the calibrated data in DATA.
        """
        if self.previous_selfcal is not None or self.input_caltable:
            return
        if is_column_in_ms(self.visfile, "CORRECTED_DATA"
                           ) and not columns_match(self.visfile, "DATA", "CORRECTED_DATA"):
            warnings.warn(
                "CORRECTED_DATA differs from DATA in {0}. gaincal solves against DATA, without the "
                "calibration in CORRECTED_DATA: split the calibrated data first "
                "(split(..., datacolumn='corrected')).".format(self.visfile)
            )

    def _require_model(self) -> None:
        """
        Protected method that raises if the imager did not write MODEL_DATA. Without it gaincal
        silently solves against a 1 Jy point source at the phase centre.
        """
        if not is_column_in_ms(self.visfile, "MODEL_DATA"):
            raise RuntimeError(
                "{0} did not write MODEL_DATA in {1}. Self-calibration solves against the model, "
                "so the imager must save it to the model column.".format(
                    type(self.imager).__name__, self.visfile
                )
            )

    def _flagged_fraction(self) -> float:
        """
        Protected method that returns the fraction of the selected data that is flagged

        Returns
        -------
        The flagged fraction, between 0 and 1
        """
        summary = flagdata(vis=self.visfile, mode="summary", field=self.field, spw=self.spw)
        return summary["flagged"] / summary["total"] if summary["total"] else 0.0

    def _report_flagged_fraction(self) -> None:
        """
        Protected method that records and prints the flagged fraction after applying an iteration's
        solutions. With applymode="calflag" the data whose solutions fail minsnr are flagged, and
        these flags add up over iterations (#27).
        """
        self._flagged_history.append(self._flagged_fraction())
        print(
            "Flagged: {0:0.2%} ({1:+0.2%} in this iteration)".format(
                self._flagged_history[-1], self._flagged_history[-1] - self._flagged_history[-2]
            )
        )

    def _has_phase_solution(self) -> bool:
        """
        Protected method that checks whether phase solutions are available to pre-apply, either
        from an input calibration table or from the chain of previous self-calibrations

        Returns
        -------
        True if phase solutions are available, False otherwise
        """
        if self.input_caltable:
            return True
        previous = self.previous_selfcal
        return previous is not None and any(
            caltable.solves_phase for caltable in previous._caltable_stack
        )

    def _require_phase_solution(self) -> None:
        """
        Protected method that raises if amplitude self-calibration would run without phase
        solutions. Phase errors decorrelate the visibilities within the solution interval and bias
        the amplitude solutions low.
        """
        if not self._has_phase_solution():
            raise ValueError(
                "{0} needs phase solutions to pre-apply: run a Phasecal first and pass it as "
                "previous_selfcal (after calling its run()), or pass an input_caltable with phase "
                "solutions.".format(type(self).__name__)
            )

    def _copy_directory_at_start(self):
        if self.visfile is not None:
            path_object = Path(self.visfile)

            current_visfile = "{0}_{2}{1}".format(
                Path.joinpath(path_object.parent, path_object.stem), path_object.suffix,
                self._calmode + "0"
            )
            # Copying dataset and overwriting if it has already been created
            if os.path.exists(current_visfile):
                shutil.rmtree(current_visfile)
            shutil.copytree(self.visfile, current_visfile)
            self.visfile = current_visfile
            self.imager.inputvis = current_visfile

    def _copy_directory_during_iterations(self, iteration):
        path_object = Path(self.visfile)

        current_visfile = "{0}_{2}{1}".format(
            Path.joinpath(path_object.parent, path_object.stem), path_object.suffix,
            self._calmode + str(iteration + 1)
        )
        # Copying dataset and overwriting if it has already been created
        if os.path.exists(current_visfile):
            shutil.rmtree(current_visfile)
        shutil.copytree(self.visfile, current_visfile)

        return current_visfile

    def _save_selfcal(self, caltable_version="", overwrite=True) -> None:
        """
        Protected function that saves the flags using CASA flag manager

        Parameters
        ----------
        caltable_version :
            Calibration table version
        overwrite :
            Whether to overwrite the .flags file or not

        Returns
        -------
        None
        """
        if overwrite:
            flagmanager(vis=self.visfile, mode='delete', versionname=caltable_version)
        flagmanager(vis=self.visfile, mode='save', versionname=caltable_version)

    def _reset_selfcal(self, caltable_version="") -> None:
        """
        Protected function that resets the flags and deletes the model column if it is present in the measurement set

        Parameters
        ----------
        caltable_version :
            Calibration table version

        Returns
        -------
        None
        """
        flagmanager(vis=self.visfile, mode='restore', versionname=caltable_version)
        clearcal(self.visfile)
        delmod(vis=self.visfile, otf=True, scr=True)

    def _restore_selfcal(self, caltable_version="") -> None:
        """
        Protected function that restores the flags of dataset to a certain version

        Parameters
        ----------
        caltable_version :
            Calibration table version

        Returns
        -------
        None
        """
        flagmanager(vis=self.visfile, mode='restore', versionname=caltable_version)
        delmod(vis=self.visfile, otf=True, scr=True)

    def _init_selfcal(self) -> None:
        """
        Protected function that initializes the input calibration tables and the PSNR history if any
        previous self-calibration object is passed to the current object

        Returns
        -------
        None
        """
        if self.previous_selfcal is not None:
            if self.previous_selfcal._caltables:
                self.input_caltable = self.previous_selfcal._caltables[-1]
            else:
                self.input_caltable = ""
            self._caltable_stack = list(self.previous_selfcal._caltable_stack)
            self._psnr_history = copy.deepcopy(self.previous_selfcal._psnr_history)
        elif self.input_caltable:
            # A user table is assumed to hold phase solutions, applied with the user spwmap
            self._caltable_stack = [CalTable(self.input_caltable, "p", tuple(self.spwmap))]

    def _preapplied_caltables(self) -> list:
        """
        Protected method that returns the tables applied on the fly while solving. By default
        self-calibration is incremental: every accepted table is pre-applied.

        Returns
        -------
        List of calibration tables in application order
        """
        return list(self._caltable_stack)

    def _applied_caltables(self, caltable: CalTable) -> list:
        """
        Protected method that returns the tables applied to the data after solving a new table

        Parameters
        ----------
        caltable :
            The newly solved calibration table

        Returns
        -------
        List of calibration tables in application order
        """
        return self._preapplied_caltables() + [caltable]

    def _accept_caltable(self, caltable: CalTable) -> None:
        """
        Protected method that records a calibration table whose iteration was accepted

        Parameters
        ----------
        caltable :
            The accepted calibration table
        """
        self._caltable_stack.append(caltable)

    def _solution_spwmap(self) -> tuple:
        """
        Protected method that returns the spectral window map needed to apply the solutions of the
        current iteration, derived from combine and the solved spectral windows. An explicit
        spwmap is respected for solutions that combine spectral windows.

        Returns
        -------
        The spectral window map as a tuple
        """
        derived = derive_spwmap(
            self.combine, get_spw_ids(self.visfile, self.spw), count_spectral_windows(self.visfile)
        )
        if derived and self.spwmap:
            if list(self.spwmap) != derived:
                warnings.warn(
                    "spwmap {0} differs from {1}, derived from combine='{2}' and spw='{3}'. "
                    "Solutions combined over spectral windows are stored under the lowest "
                    "selected spectral window.".format(
                        list(self.spwmap), derived, self.combine, self.spw
                    )
                )
            return tuple(self.spwmap)
        return tuple(derived)

    def _solve(self, caltable: str, solint: str, **kwargs) -> CalTable:
        """
        Protected method that solves the gains of the current iteration with gaincal, pre-applying
        the tables returned by _preapplied_caltables

        Parameters
        ----------
        caltable :
            Path of the output calibration table
        solint :
            Solution interval
        kwargs :
            Extra gaincal arguments of the calibration mode, e.g. solnorm

        Returns
        -------
        The solved calibration table
        """
        gaincal(
            vis=self.visfile,
            caltable=caltable,
            field=self.field,
            spw=self.spw,
            uvrange=self.uvrange,
            gaintype=self.gaintype,
            refant=self.refant,
            calmode=self._calmode,
            combine=self.combine,
            solint=solint,
            minsnr=self.minsnr,
            minblperant=self.minblperant,
            **build_gaintable_arguments(self._preapplied_caltables(), self.interp),
            **kwargs
        )
        return CalTable(caltable, self._calmode, self._solution_spwmap())

    def _apply(self, caltable: CalTable) -> None:
        """
        Protected method that applies the tables returned by _applied_caltables to the data

        Parameters
        ----------
        caltable :
            The newly solved calibration table
        """
        print("Applying calibration tables to {0} file".format(self.visfile))
        applycal(
            vis=self.visfile,
            field=self.field,
            spw=self.spw,
            gainfield='',
            calwt=False,
            flagbackup=False,
            applymode=self.applymode,
            **build_gaintable_arguments(self._applied_caltables(caltable), self.interp)
        )

    def _init_run(self, image_name_string: str = "") -> None:
        """
        Protected function that runs the imager at the beginning of the self-calibration run in order to initializes
        the model column

        Parameters
        ----------
        image_name_string :
            The string of the resulting CASA image name
        """
        self._warn_if_corrected_differs_from_data()
        self._flagged_history.append(self._flagged_fraction())
        print("Flagged before self-calibration: {0:0.2%}".format(self._flagged_history[-1]))
        if not self._ismodel_in_dataset() or self.previous_selfcal is None:
            if self.previous_selfcal is None:
                # A virtual model left in the dataset (e.g. by setjy or tclean
                # savemodel="virtual") would compete with the model column the imager writes
                delmod(vis=self.visfile, otf=True, scr=False)
            imagename = self._image_name + image_name_string
            self.imager.run(imagename)
            self._require_model()
            print("Original: - PSNR: {0:0.3f}".format(self.imager.psnr))
            print("Peak: {0:0.3f} mJy/beam".format(self.imager.peak * 1000.0))
            print("Noise: {0:0.3f} mJy/beam".format(self.imager.stdv * 1000.0))
            self._psnr_history.append(self.imager.psnr)

    def _run_imager(self, current_iteration: int = 0) -> None:
        """
        Protected method that runs the imager at a certain self-calibration iteration

        Parameters
        ----------
        iter :
            Iteration number during the loop
        """
        self._report_flagged_fraction()

        imagename = self._image_name + '_' + self._calmode + str(current_iteration)

        self.imager.run(imagename)
        self._require_model()

        self._psnr_history.append(self.imager.psnr)

        print(
            "Solint: {0} - PSNR: {1:0.3f}".format(
                self.solint[current_iteration], self._psnr_history[-1]
            )
        )
        print("Peak: {0:0.3f} mJy/beam".format(self.imager.peak * 1000.0))
        print("Noise: {0:0.3f} mJy/beam".format(self.imager.stdv * 1000.0))

    def _finish_selfcal_iteration(self, current_iteration: int = 0) -> bool:
        """
        Protected method that finishes self-calibration iterations. If the PSNR of the current iteration improves then
        a new dataset is created and the measurement set file name is changed. Otherwise the flags are restored to the
        last version and the PSNR history and last calibration table are popped from the lists. The measurement set
        name is changed to the last (the one that had better PSNR).

        Parameters
        ----------
        current_iteration :
            Iteration number during the self-calibration loop

        Returns
        -------

        """
        if self.restore_psnr:
            if len(self._psnr_history) > 1:

                print(
                    "Old PSNR {0:0.3f} vs last PSNR {1:0.3f}".format(
                        self._psnr_history[-2], self._psnr_history[-1]
                    )
                )

                if self._psnr_history[-1] <= self._psnr_history[-2]:

                    print(
                        "PSNR decreasing or equal in this solution interval - restoring to last MS and exiting loop..."
                    )
                    self._restore_selfcal(caltable_version=self._caltables_versions[-1])
                    self._psnr_history.pop()
                    self._caltables.pop()
                    # Restoring to last MS
                    self.visfile = self._psnr_visfile_backup
                    self.imager.inputvis = self._psnr_visfile_backup
                    return True
                else:
                    print(
                        "PSNR improved on iteration {0} - Copying measurement set files...".
                        format(current_iteration)
                    )

                    if current_iteration + 1 < self._loops:
                        current_visfile = self._copy_directory_during_iterations(current_iteration)

                        # Saving old visfile name
                        self._psnr_visfile_backup = self.visfile
                        # Changing visfile attribute to new current_visfile for selfcal and imager
                        self.visfile = current_visfile
                        self.imager.inputvis = current_visfile
                    return False
            else:
                return False
        else:
            return False

    def _flag_dataset(
        self, datacolumn=None, mode="rflag", timedevscale=3.0, freqdevscale=3.0
    ) -> None:
        """
        Protected method that flag the dataset on each iteration. This functions aims to flag residual outliers.

        NOTE 1: RESIDUAL = CORRECTED - MODEL
        RESIDUAL_DATA = DATA - MODEL

        NOTE 2: When datacolumn is WEIGHT, the task will internally use WEIGHT_SPECTRUM.
        If WEIGHT_SPECTRUM does not exist, it will create one on-the-fly based on the values of WEIGHT.

        Parameters
        ----------
        datacolumn :
            The datacolumn to use in order to flag. (See description).
        mode :
            Flagging mode, one of FLAG_MODES
        timedevscale :
            rflag only. For time analysis, flag a point if local RMS around it is larger than timedevscale $x$ timedev.
        freqdevscale :
            rflag only. For spectral analysis, flag a point if local rms around it is larger than freqdevscale $x$ freqdev.
        """

        if datacolumn is None:
            if is_column_in_ms(self.visfile, "CORRECTED_DATA"):
                datacolumn = "residual"
            else:
                datacolumn = "residual_data"

        # Only the parameters of the chosen mode: flagdata ignores the others (#27)
        if mode == "rflag":
            mode_parameters = dict(
                timedevscale=timedevscale, freqdevscale=freqdevscale, spectralmax=500
            )
        else:
            mode_parameters = dict(
                timecutoff=5.0, freqcutoff=5.0, freqfit='line', flagdimension='freq'
            )

        print("Flagging {0} data column using {1}".format(datacolumn, mode))

        # The flags of each iteration are saved by _save_selfcal, so no extra flag backup is made
        flagdata(
            vis=self.visfile,
            mode=mode,
            datacolumn=datacolumn,
            field=self.field,
            spw=self.spw,
            extendflags=False,
            ntime="scan",
            action='apply',
            flagbackup=False,
            **mode_parameters
        )

    def _ismodel_in_dataset(self) -> bool:
        """
        Protected function that checks if the MODEL_DATA column exists in the current measurement set file

        Returns
        -------
        True if the MODEL_DATA column exists, False otherwise
        """
        return is_column_in_ms(self.visfile, "MODEL_DATA")

    def _set_attributes_from_dicts(self, current_iteration: int = 0) -> None:
        """
        Protected method that assigns iteration attributes from dictionaries to imager and this self-calibration
        instance.

        Parameters
        ----------
        current_iteration :
            Current iteration in the self-calibration loop
        """
        changes = ((self.imager, self.varchange_imager), (self, self.varchange_selfcal))
        for target, varchange in changes:
            for key, values in (varchange or {}).items():
                # A selection change applies to the solve, the apply and the imager alike
                for changed in (self, self.imager) if key in SELECTION_PARAMETERS else (target, ):
                    setattr(changed, key, values[current_iteration])

    def _plot_selfcal(
        self,
        caltable,
        xaxis="",
        yaxis="",
        iteration="",
        timerange="",
        antenna="",
        subplot=[1, 1],
        plotrange=[],
        want_plot=False,
        **kwargs
    ):
        figfile_name = caltable + ".png"
        if want_plot:
            print("Plot selfcal is trying to plot, but the function is not part of CASA 6 yet")
            # au.plotms(vis=caltable, xaxis=xaxis, yaxis=yaxis, iteration=iteration, gridrows=subplot[0],
            # gridcols=subplot[1], antenna=antenna, timerange=timerange, plotrange=plotrange, plotfile=figfile_name,
            # overwrite=True, showgui=True)
        else:
            print("Plot selfcal is trying to plot, but the function is not part of CASA 6 yet")
            # au.plotms(vis=caltable, xaxis=xaxis, yaxis=yaxis, iteraxis=iteration, gridrows=subplot[0],
            # gridcols=subplot[1], antenna=antenna, timerange=timerange, plotrange=plotrange, plotfile=figfile_name,
            # overwrite=True, showgui=False)

    def selfcal_output(self, overwrite=False, _statwt=False, min_samp=8) -> str:
        """
        Public function that creates a new measurement set only taking the corrected column.
        If _statwt is True then applies the statwt function and creates a .statwt measurement

        Parameters
        ----------
        overwrite :
            Whether to overwrite the measurement set files
        _statwt :
            Whether to create a new measurement applying the statwt function
        min_samp :
            Minimum number of unflagged visibilities for estimating the scatter if _statwt is True
        Returns
        -------
        Name of the self-calibrated measurement set file
        """
        output_vis = self.visfile + '.selfcal'
        data_column = ""
        if overwrite:
            if os.path.exists(output_vis):
                shutil.rmtree(output_vis)
        if is_column_in_ms(self.visfile, "CORRECTED_DATA"):
            data_column = "corrected"
        else:
            # Raise warning
            warnings.warn(
                "Corrected data column is not present, data column will be extracted instead."
            )
            data_column = "data"
        split(vis=self.visfile, outputvis=output_vis, datacolumn=data_column)
        if _statwt:
            statwt_path = output_vis + '.statwt'
            if os.path.exists(statwt_path):
                shutil.rmtree(statwt_path)
            shutil.copytree(output_vis, statwt_path)
            statwt(vis=statwt_path, datacolumn="data", minsamp=min_samp)
        return output_vis

    def _uvsubtract(self):
        uvsub(vis=self.visfile, reverse=False)

    def _uvadd(self):
        uvsub(vis=self.visfile, reverse=True)

    @abstractmethod
    def run(self):
        """
        Abstract method that runs the self-calibration
        """
        pass
