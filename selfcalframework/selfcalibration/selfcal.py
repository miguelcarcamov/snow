import os
import shutil
import sys
from abc import ABCMeta, abstractmethod
from datetime import datetime
from pathlib import Path

from casatasks import (
    applycal, clearcal, delmod, flagdata, flagmanager, gaincal, rmtables, split, statwt, uvsub
)
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
        self._psnr_file_backup = ""
        self._calmode = ""

        if output_caltables is None:
            self.output_caltables = self.imager.output

        if self.imager is None:
            print("Error, Imager Object is Nonetype")
            raise ValueError("Error, self-calibration objects cannot run without an imager object")

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

        if self.subtract_source == True:
            if self.Imager.getPhaseCenter() != "":
                raise ValueError(
                    "Error, phase center needs to be set if a source is going to be subtracted"
                )

    @staticmethod
    def read_last_line_file_backup(file_name: str = ""):
        first_line = ""
        with open(file_name, 'rb') as f:
            file.seek(-2, os.SEEK_END)
            while f.read(1) != b'\n':
                f.seek(-2, os.SEEK_CUR)
            return float(f.readline().decode())

    def write_file_backup(self):
        with open(self.psnr_file_backup, 'a') as f:
            f.write("{0:0.5f}\n".format(self.Imager.getPSNR()))

    def delete_last_lines(self):
        with open(self.psnr_file_backup, 'r+') as f:
            lines = f.readlines()
            fp.seek(0)
            fp.truncate()
            fp.writelines(lines[0:1])

    def save_selfcal(self, caltable_version="", overwrite=True):
        if overwrite:
            flagmanager(vis=self.visfile, mode='delete', versionname=caltable_version)
        flagmanager(vis=self.visfile, mode='save', versionname=caltable_version)

    def reset_selfcal(self, caltable_version=""):
        flagmanager(vis=self.visfile, mode='restore', versionname=caltable_version)
        clearcal(self.visfile)
        delmod(vis=self.visfile, otf=True, scr=True)

    def restore_selfcal(self, caltable_version=""):
        flagmanager(vis=self.visfile, mode='restore', versionname=caltable_version)
        delmod(vis=self.visfile, otf=True, scr=True)

    def flag_dataset(self, datacolumn="RESIDUAL", mode="rflag", timedevscale=3.0, freqdevscale=3.0):
        # NOTE1: RESIDUAL = CORRECTED - MODEL
        # RESIDUAL_DATA = DATA - MODEL
        # NOTE2: When datacolumn is WEIGHT, the task will
        # internally use WEIGHT_SPECTRUM.
        # If WEIGHT_SPECTRUM does not exist, it will create one on-the-fly based on the values of WEIGHT.
        flagdata(
            vis=self.visfile,
            mode=mode,
            datacolumn=datacolumn,
            field='',
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

    def ismodel_in_dataset(self):
        tb.open(tablename=self.visfile)
        columns = tb.colnames()
        tb.close()
        if "MODEL_DATA" in columns:
            return True
        else:
            return False

    def set_attributes_from_dicts(self, iteration=0):
        if self.varchange_imager is not None:
            for key in self.varchange_imager.keys():
                setattr(self.Imager, key, self.varchange_imager[key][iteration])

        if self.varchange_selfcal is not None:
            for key in self.varchange_selfcal.keys():
                setattr(self, key, self.varchange_selfcal[key][iteration])

    def plot_selfcal(
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
        outputvis = self.visfile + '.selfcal'

        if overwrite:
            os.system('rm -rf ' + outputvis)
        split(vis=self.visfile, outputvis=outputvis, datacolumn='corrected')
        if _statwt:
            statwt_path = outputvis + '.statwt'
            if os.path.exists(statwt_path):
                shutil.rmtree(statwt_path)
            shutil.copytree(outputvis, statwt_path)
            statwt(vis=statwt_path, datacolumn="data")
        return outputvis

    def uvsubtract(self):
        uvsub(vis=self.visfile, reverse=False)

    def uvsubtract(vis=""):
        uvsub(vis=vis, reverse=False)

    def uvadd(self):
        uvsub(vis=self.visfile, reverse=True)

    def uvadd(vis=""):
        uvsub(vis=vis, reverse=True)

    @abstractmethod
    def run(self):
        return
