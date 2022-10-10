from __future__ import annotations

import copy
import os
import shutil
from abc import ABCMeta, abstractmethod
from pathlib import Path
from dataclasses import dataclass

from casatasks import (clearcal, delmod, flagdata, flagmanager, split, statwt, uvsub)
from casatools import table

from ..imaging.imager import Imager

tb = table()


@dataclass(init=False, repr=True)
class Selfcal(metaclass=ABCMeta):

    def __init__(
        self,
        visfile: str = "",
        imager: Imager = None,
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
            Imager instance to run
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
            Calibration mode - ””=”calflag”, ”calflagstrict”, ”trial”, ”flagonly”, ”flagonlystrict”, or ”calonly”
        flag_mode :
            Flag mode operation. e.g : "manual", "clip", "quack", "shadow", "elevation", "tfcrop", "rflag"
        combine :
            Data axes to combine for solving
        flag_dataset :
            Whether to flag Fourier residuals outliers
        restore_psnr :
            Restores the dataset if the peak signal-to-noise ratio decreases
        subtract_source :
            Subtract source model if needed
        """
        # Public variables
        self.visfile = visfile
        self.imager = imager
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
        self._caltables_versions = []
        self._psnr_history = []
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

    def _copy_directory_at_start(self):
        if self.visfile is not None:
            path_object = Path(self.visfile)

            current_visfile = "{0}_{2}{1}".format(
                Path.joinpath(path_object.parent, path_object.stem), path_object.suffix,
                self._calmode
            )
            # Copying dataset and overwriting if it has already been created
            if os.path.exists(current_visfile):
                shutil.rmtree(current_visfile)
            shutil.copytree(self.visfile, current_visfile)
            self.visfile = current_visfile

    def _copy_directory_during_iterations(self, iteration):
        path_object = Path(self.visfile)

        current_visfile = "{0}_{2}{1}".format(
            Path.joinpath(path_object.parent, path_object.stem), path_object.suffix, str(iteration)
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
            self._psnr_history = copy.deepcopy(self.previous_selfcal._psnr_history)

    def _init_run(self, image_name_string: str = "") -> None:
        """
        Protected function that runs the imager at the beginning of the self-calibration run in order to initializes
        the model column

        Parameters
        ----------
        image_name_string :
            The string of the resulting CASA image name
        """
        if not self._ismodel_in_dataset() or self.previous_selfcal is None:
            imagename = self._image_name + image_name_string
            self.imager.run(imagename)
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
        imagename = self._image_name + '_' + self._calmode + str(current_iteration)

        self.imager.run(imagename)

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
                    return True
                else:
                    print(
                        "PSNR improved on iteration {0} - Copying measurement set files...".
                        format(current_iteration)
                    )

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
        self, datacolumn="RESIDUAL", mode="rflag", timedevscale=3.0, freqdevscale=3.0
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
            Flagging mode
        timedevscale :
            For time analysis, flag a point if local RMS around it is larger than timedevscale $x$ timedev.
        freqdevscale :
            For spectral analysis, flag a point if local rms around it is larger than freqdevscale $x$ freqdev.
        """

        flagdata(
            vis=self.visfile,
            mode=mode,
            datacolumn=datacolumn,
            field=self.imager.field,
            timecutoff=5.0,
            freqcutoff=5.0,
            timefit='line',
            freqfit='line',
            flagdimension='freqtime',
            extendflags=False,
            timedevscale=timedevscale,
            freqdevscale=freqdevscale,
            spectralmax=500,
            extendpols=False,
            growaround=False,
            flagneartime=False,
            flagnearfreq=False,
            action='apply',
            flagbackup=True,
            overwrite=True,
            writeflags=True
        )

    def _ismodel_in_dataset(self) -> bool:
        """
        Protected function that checks if the MODEL_DATA column exists in the current measurement set file

        Returns
        -------
        None
        """
        tb.open(tablename=self.visfile)
        columns = tb.colnames()
        tb.close()
        if "MODEL_DATA" in columns:
            return True
        else:
            return False

    def _set_attributes_from_dicts(self, current_iteration: int = 0) -> None:
        """
        Protected method that assigns iteration attributes from dictionaries to imager and this self-calibration
        instance.

        Parameters
        ----------
        current_iteration :
            Current iteration in the self-calibration loop
        """
        if self.varchange_imager is not None:
            for key in self.varchange_imager.keys():
                setattr(self.imager, key, self.varchange_imager[key][current_iteration])

        if self.varchange_selfcal is not None:
            for key in self.varchange_selfcal.keys():
                setattr(self, key, self.varchange_selfcal[key][current_iteration])

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

    def selfcal_output(self, overwrite=False, _statwt=False) -> str:
        """
        Public function that creates a new measurement set only taking the corrected column.
        If _statwt is True then applies the statwt function and creates a .statwt measurement

        Parameters
        ----------
        overwrite :
            Whether to overwrite the measurement set files
        _statwt :
            Whether to create a new measurement applying the statwt function

        Returns
        -------
        Name of the self-calibrated measurement set file
        """
        output_vis = self.visfile + '.selfcal'

        if overwrite:
            os.system('rm -rf ' + output_vis)
        split(vis=self.visfile, outputvis=output_vis, datacolumn='corrected')
        if _statwt:
            statwt_path = output_vis + '.statwt'
            if os.path.exists(statwt_path):
                shutil.rmtree(statwt_path)
            shutil.copytree(output_vis, statwt_path)
            statwt(vis=statwt_path, datacolumn="data")
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
