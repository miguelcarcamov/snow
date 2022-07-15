from __future__ import annotations

import copy
import os
import shutil
from abc import ABCMeta, abstractmethod
from pathlib import Path

from casatasks import (clearcal, delmod, flagdata, flagmanager, split, statwt, uvsub)
from casatools import table

from ..imaging.imager import Imager

tb = table()


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
        self._psnr_visfile_backup = ""
        self._calmode = ""
        self._loops = 0

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

    def _save_selfcal(self, caltable_version="", overwrite=True):
        if overwrite:
            flagmanager(vis=self.visfile, mode='delete', versionname=caltable_version)
        flagmanager(vis=self.visfile, mode='save', versionname=caltable_version)

    def _reset_selfcal(self, caltable_version=""):
        flagmanager(vis=self.visfile, mode='restore', versionname=caltable_version)
        clearcal(self.visfile)
        delmod(vis=self.visfile, otf=True, scr=True)

    def _restore_selfcal(self, caltable_version=""):
        flagmanager(vis=self.visfile, mode='restore', versionname=caltable_version)
        delmod(vis=self.visfile, otf=True, scr=True)

    def _init_selfcal(self):
        if self.previous_selfcal is not None:
            if self.previous_selfcal._caltables:
                self.input_caltable = self.previous_selfcal._caltables[-1]
            else:
                self.input_caltable = ""
            self._psnr_history = copy.deepcopy(self.previous_selfcal._psnr_history)

    def _init_run(self, image_name_string: str = ""):
        if not self._ismodel_in_dataset() or self.previous_selfcal is None:
            imagename = self._image_name + image_name_string
            self.imager.run(imagename)
            print("Original: - PSNR: {0:0.3f}".format(self.imager.psnr))
            print("Peak: {0:0.3f} mJy/beam".format(self.imager.peak * 1000.0))
            print("Noise: {0:0.3f} mJy/beam".format(self.imager.stdv * 1000.0))
            self._psnr_history.append(self.imager.psnr)

    def _run_imager(self, iter: int = 0):
        imagename = self._image_name + '_' + self._calmode + str(iter)

        self.imager.run(imagename)

        self._psnr_history.append(self.imager.psnr)

        print("Solint: {0} - PSNR: {1:0.3f}".format(self.solint[iter], self._psnr_history[-1]))
        print("Peak: {0:0.3f} mJy/beam".format(self.imager.peak * 1000.0))
        print("Noise: {0:0.3f} mJy/beam".format(self.imager.stdv * 1000.0))

    def _finish_selfcal_iteration(self, iter: int = 0):

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
                        format(iter)
                    )

                    path_object = Path(self.visfile)

                    current_visfile = "{0}_{2}{1}".format(
                        Path.joinpath(path_object.parent, path_object.stem), path_object.suffix,
                        self._calmode + str(iter)
                    )
                    # Copying dataset and overwriting if it has already been created
                    if os.path.exists(current_visfile):
                        shutil.rmtree(current_visfile)
                    shutil.copytree(self.visfile, current_visfile)
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
    ):
        # NOTE1: RESIDUAL = CORRECTED - MODEL
        # RESIDUAL_DATA = DATA - MODEL
        # NOTE2: When datacolumn is WEIGHT, the task will
        # internally use WEIGHT_SPECTRUM.
        # If WEIGHT_SPECTRUM does not exist, it will create one on-the-fly based on the values of WEIGHT.
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

    def _ismodel_in_dataset(self):
        tb.open(tablename=self.visfile)
        columns = tb.colnames()
        tb.close()
        if "MODEL_DATA" in columns:
            return True
        else:
            return False

    def _set_attributes_from_dicts(self, iteration=0):
        if self.varchange_imager is not None:
            for key in self.varchange_imager.keys():
                setattr(self.imager, key, self.varchange_imager[key][iteration])

        if self.varchange_selfcal is not None:
            for key in self.varchange_selfcal.keys():
                setattr(self, key, self.varchange_selfcal[key][iteration])

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

    def selfcal_output(self, overwrite=False, _statwt=False):
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

    def _uvsubtract(self):
        uvsub(vis=self.visfile, reverse=False)

    def _uvadd(self):
        uvsub(vis=self.visfile, reverse=True)

    @abstractmethod
    def run(self):
        return
