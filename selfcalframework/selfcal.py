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

tb = table()


class Selfcal(metaclass=ABCMeta):

    def __init__(
        self,
        visfile="",
        Imager=None,
        refant="",
        spwmap=[],
        minblperant=4,
        want_plot=True,
        interp='linear',
        gaintype='T',
        uvrange="",
        solint=[],
        varchange_imager=None,
        varchange_selfcal=None,
        output_caltables="",
        minsnr=3.0,
        applymode="calflag",
        flag_mode="rflag",
        combine="",
        flag_dataset_bool=False,
        restore_PSNR=False,
        subtract_source=False
    ):
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.caltables = []
        self.caltables_versions = []
        self.psnr_history = []
        self.output_caltables = self.Imager.output
        self.psnr_file_backup = ""

        if self.Imager is None:
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

    def getVisfile(self):
        return self.getvisfile

    def getImager(self):
        return self.Imager

    def getRefant(self):
        return self.refant

    def getSpwmap(self):
        return self.spwmap

    def getCaltables(self):
        return self.caltables

    def getCaltablesVersions(self):
        return self.caltables_versions

    def getPSNRHistory(self):
        return self.psnr_history

    def getSubtractSource(self):
        return self.subtract_source

    def read_first_line_file_backup(self):
        first_line = ""
        with open(self.psnr_file_backup, 'r') as f:
            first_line = f.readline().rstrip()
        return float(first_line)

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


class Ampcal(Selfcal):

    def __init__(self, selfcal_object=None, input_caltable="", **kwargs):
        super(Ampcal, self).__init__(**kwargs)
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])

        self.calmode = 'a'
        self.loops = len(self.solint)
        self.imagename = self.Imager.getOutputPath()
        self.psnr_file_backup = self.Imager.output + "psnr_amp.txt"

        if self.selfcal_object is None and self.input_caltable == "":
            print("Error, Self-cal object is Nonetype and input_caltable is an empty string")
            raise ValueError(
                "Error, Amplitude self-cal objects cannot run without an phase-cal object or input caltable"
            )
        else:
            if self.selfcal_object:
                self.input_caltable = self.selfcal_object.getCaltables()[-1]
            elif self.input_caltable != "":
                if not os.path.exists(self.input_caltable):
                    raise FileNotFoundError(
                        "The caltable " + self.input_caltable + " needs to be created"
                    )
            else:
                print("Error, Ampcal needs a non-empty list of caltables")
                raise ValueError(
                    "Error, Amplitude self-cal objects cannot run with an empty caltable list"
                )

    def run(self):
        caltable = ""
        if not self.ismodel_in_dataset():
            imagename = self.imagename + "before_apcal"
            self.Imager.run(imagename)
            print("Before amplitude self-cal: - PSNR: " + str(self.Imager.getPSNR()))
            print("Noise: " + str(self.Imager.getSTDV() * 1000.0) + " mJy/beam")
            self.write_file_backup()
            self.psnr_history.append(self.Imager.getPSNR())

            path_object = Path(self.visfile)
            current_visfile = "{0}_{2}{1}".format(
                Path.joinpath(path_object.parent, path_object.stem), path_object.suffix,
                "before_apcal"
            )
            if os.path.exists(current_visfile):
                shutil.rmtree(current_visfile)
            shutil.copytree(self.visfile, current_visfile)
        else:
            psnr_file = self.read_first_line_file_backup()
            self.psnr_history.append(psnr_file)
            self.delete_last_lines()

        for i in range(0, self.loops):
            caltable = self.output_caltables + 'ampcal_' + str(i)
            self.caltables.append(caltable)
            rmtables(caltable)

            self.set_attributes_from_dicts(i)

            gaincal(
                vis=self.visfile,
                field=self.Imager.getField(),
                caltable=caltable,
                spw=self.Imager.getSpw(),
                uvrange=self.uvrange,
                gaintype=self.gaintype,
                refant=self.refant,
                calmode=self.calmode,
                combine=self.combine,
                solint=self.solint[i],
                minsnr=self.minsnr,
                minblperant=self.minblperant,
                gaintable=self.input_caltable,
                spwmap=self.spwmap,
                solnorm=True
            )

            self.plot_selfcal(
                caltable,
                xaxis="time",
                yaxis="amp",
                iteration="antenna",
                subplot=[4, 2],
                plotrange=[0, 0, 0.2, 1.8],
                want_plot=self.want_plot
            )

            versionname = 'before_ampcal_' + str(i)
            self.save_selfcal(caltable_version=versionname, overwrite=True)
            self.caltables_versions.append(versionname)
            applycal(
                vis=self.visfile,
                spw=self.Imager.getSpw(),
                spwmap=[self.spwmap, self.spwmap],
                field=self.Imager.getField(),
                gaintable=[self.input_caltable, caltable],
                gainfield='',
                calwt=False,
                flagbackup=False,
                interp=self.interp,
                applymode=self.applymode
            )

            if self.flag_dataset_bool:
                self.flag_dataset(mode=self.flag_mode)

            imagename = self.imagename + '_a' + str(i)

            self.Imager.run(imagename)

            self.psnr_history.append(self.Imager.getPSNR())

            print("Solint: " + str(self.solint[i]) + " - PSNR: " + str(self.psnr_history[i]))
            print("Noise: " + str(self.Imager.getSTDV() * 1000.0) + " mJy/beam")

            path_object = Path(self.visfile)
            current_visfile = "{0}_{2}{1}".format(
                Path.joinpath(path_object.parent, path_object.stem), path_object.suffix,
                "amp" + str(i)
            )
            if self.restore_PSNR:
                if i > 0:
                    if self.psnr_history[-1] <= self.psnr_history[-2]:
                        self.restore_selfcal(caltable_version=self.caltables_versions[-1])
                        self.psnr_history.pop()
                        self.caltables.pop()
                        print(
                            "PSNR decreasing in this solution interval - restoring to last MS and exiting loop"
                        )
                        break
                    else:
                        print(
                            "PSNR improved on iteration {0} - Copying measurement set files...".
                            format(i)
                        )
                        if os.path.exists(current_visfile):
                            shutil.rmtree(current_visfile)
                        shutil.copytree(self.visfile, current_visfile)
                        self.visfile = current_visfile
                        self.Imager.inputvis = current_visfile
                        self.write_file_backup()

                elif self.selfcal_object:
                    if self.psnr_history[i] <= self.selfcal_object.getPSNRHistory()[-1]:
                        self.restore_selfcal(
                            caltable_version=self.selfcal_object.getCaltablesVersions()[-1]
                        )
                        self.psnr_history.pop()
                        self.caltables.pop()
                        self.caltables = self.selfcal_object.getCaltables()
                        self.psnr_history = self.selfcal_object.getPSNRHistory()
                        self.caltables_versions = self.selfcal_object.getCaltablesVersions()
                        print(
                            "PSNR decreasing or equal in this solution interval - restoring to last MS and "
                            "exiting loop"
                        )
                        break
                    else:
                        print(
                            "PSNR improved on iteration {0} - Copying measurement set files...".
                            format(i)
                        )
                        if os.path.exists(current_visfile):
                            shutil.rmtree(current_visfile)
                        shutil.copytree(self.visfile, current_visfile)
                        self.visfile = current_visfile
                        self.Imager.inputvis = current_visfile
                        self.write_file_backup()


class Phasecal(Selfcal):

    def __init__(self, input_caltable="", **kwargs):
        super(Phasecal, self).__init__(**kwargs)
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])

        self.calmode = 'p'
        self.loops = len(self.solint)
        self.imagename = self.Imager.getOutputPath()
        self.psnr_file_backup = self.Imager.output + "psnr_ph.txt"

    def run(self):
        caltable = "before_selfcal"
        self.save_selfcal(caltable_version=caltable, overwrite=True)
        self.caltables_versions.append(caltable)
        if not self.ismodel_in_dataset():
            imagename = self.imagename + '_original'
            self.Imager.run(imagename)
            print("Original: - PSNR: " + str(self.Imager.getPSNR()))
            print("Noise: " + str(self.Imager.getSTDV() * 1000.0) + " mJy/beam")
            self.write_file_backup()
            self.psnr_history.append(self.Imager.getPSNR())
        else:
            psnr_file = self.read_first_line_file_backup()
            self.psnr_history.append(psnr_file)
            self.delete_last_lines()

            path_object = Path(self.visfile)
            current_visfile = "{0}_{2}{1}".format(
                Path.joinpath(path_object.parent, path_object.stem), path_object.suffix, "original"
            )
            if os.path.exists(current_visfile):
                shutil.rmtree(current_visfile)
            shutil.copytree(self.visfile, current_visfile)

        for i in range(0, self.loops):
            caltable = 'pcal' + str(i)
            self.caltables.append(caltable)
            rmtables(caltable)

            self.set_attributes_from_dicts(i)

            gaincal(
                vis=self.visfile,
                caltable=caltable,
                field=self.Imager.getField(),
                spw=self.Imager.getSpw(),
                uvrange=self.uvrange,
                gaintype=self.gaintype,
                refant=self.refant,
                calmode=self.calmode,
                combine=self.combine,
                solint=self.solint[i],
                minsnr=self.minsnr,
                spwmap=self.spwmap,
                minblperant=self.minblperant
            )

            self.plot_selfcal(
                caltable,
                xaxis="time",
                yaxis="phase",
                iteration="antenna",
                subplot=[4, 2],
                plotrange=[0, 0, -180, 180],
                want_plot=self.want_plot
            )

            versionname = 'before_phasecal_' + str(i)
            self.save_selfcal(caltable_version=versionname, overwrite=True)
            self.caltables_versions.append(versionname)

            applycal(
                vis=self.visfile,
                field=self.Imager.getField(),
                spw=self.Imager.getSpw(),
                spwmap=self.spwmap,
                gaintable=[caltable],
                gainfield='',
                calwt=False,
                flagbackup=False,
                interp=self.interp,
                applymode=self.applymode
            )

            if self.flag_dataset_bool:
                self.flag_dataset(mode=self.flag_mode)

            imagename = self.imagename + '_ph' + str(i)

            self.Imager.run(imagename)

            self.psnr_history.append(self.Imager.getPSNR())

            print("Solint: " + str(self.solint[i]) + " - PSNR: " + str(self.psnr_history[-1]))
            print("Noise: " + str(self.Imager.getSTDV() * 1000.0) + " mJy/beam")
            path_object = Path(self.visfile)
            current_visfile = "{0}_{2}{1}".format(
                Path.joinpath(path_object.parent, path_object.stem), path_object.suffix,
                "ph" + str(i)
            )
            if os.path.exists(current_visfile):
                shutil.rmtree(current_visfile)
            shutil.copytree(self.visfile, current_visfile)
            if self.restore_PSNR:
                if len(self.psnr_history) > 1:
                    if self.psnr_history[-1] <= self.psnr_history[-2]:
                        self.restore_selfcal(caltable_version=self.caltables_versions[i])
                        self.psnr_history.pop()
                        self.caltables.pop()
                        print(
                            "PSNR decreasing or equal in this solution interval - restoring to last MS and exiting loop"
                        )
                        break
                    else:
                        print(
                            "PSNR improved on iteration {0} - Copying measurement set files...".
                            format(i)
                        )
                        if os.path.exists(current_visfile):
                            shutil.rmtree(current_visfile)
                        shutil.copytree(self.visfile, current_visfile)
                        self.visfile = current_visfile
                        self.Imager.inputvis = current_visfile
                        self.write_file_backup()


class AmpPhasecal(Selfcal):

    def __init__(
        self, selfcal_object=None, input_caltable="", incremental=False, solnorm=True, **kwargs
    ):
        super(AmpPhasecal, self).__init__(**kwargs)
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])

        self.calmode = 'ap'
        self.loops = len(self.solint)
        self.imagename = self.Imager.getOutputPath()
        self.psnr_file_backup = self.Imager.output + "psnr_ap.txt"

        if self.selfcal_object is None and self.input_caltable == "":
            raise ValueError(
                "Error, Self-cal object is Nonetype and input_caltable is an empty string"
            )
        else:
            if self.selfcal_object:
                self.input_caltable = self.selfcal_object.getCaltables()[-1]
            elif self.input_caltable != "":
                if not os.path.exists(self.input_caltable):
                    raise FileNotFoundError(
                        "The caltable " + self.input_caltable + " needs to be created"
                    )
            else:
                print("Error, Ampcal needs a non-empty list of caltables")
                sys.exit("Error, Amplitude self-cal objects cannot run with an empty caltable list")

    def run(self):
        caltable = ""
        if not self.ismodel_in_dataset():
            imagename = self.imagename + "before_apcal"
            self.Imager.run(imagename)
            print("Before amplitude-phase self-cal: - PSNR: " + str(self.Imager.getPSNR()))
            print("Noise: " + str(self.Imager.getSTDV() * 1000.0) + " mJy/beam")
            self.write_file_backup()
            self.psnr_history.append(self.Imager.getPSNR())
        else:
            psnr_file = self.read_first_line_file_backup()
            self.psnr_history.append(psnr_file)
            self.delete_last_lines()

        for i in range(0, self.loops):
            caltable = 'apcal_' + str(i)
            self.caltables.append(caltable)
            rmtables(caltable)

            self.set_attributes_from_dicts(i)

            if self.incremental:
                gaincal(
                    vis=self.visfile,
                    field=self.Imager.getField(),
                    caltable=caltable,
                    spw=self.Imager.getSpw(),
                    uvrange=self.uvrange,
                    gaintype=self.gaintype,
                    refant=self.refant,
                    calmode=self.calmode,
                    combine=self.combine,
                    solint=self.solint[i],
                    minsnr=self.minsnr,
                    minblperant=self.minblperant,
                    gaintable=self.input_caltable,
                    spwmap=self.spwmap,
                    solnorm=self.solnorm
                )
            else:
                gaincal(
                    vis=self.visfile,
                    field=self.Imager.getField(),
                    caltable=caltable,
                    spw=self.Imager.getSpw(),
                    uvrange=self.uvrange,
                    gaintype=self.gaintype,
                    refant=self.refant,
                    calmode=self.calmode,
                    combine=self.combine,
                    solint=self.solint[i],
                    minsnr=self.minsnr,
                    minblperant=self.minblperant,
                    spwmap=self.spwmap,
                    solnorm=self.solnorm
                )

            self.plot_selfcal(
                caltable,
                xaxis="time",
                yaxis="amp",
                iteration="antenna",
                subplot=[4, 2],
                plotrange=[0, 0, 0.2, 1.8],
                want_plot=self.want_plot
            )

            versionname = 'before_apcal_' + str(i)
            self.save_selfcal(caltable_version=versionname, overwrite=True)
            self.caltables_versions.append(versionname)

            if self.incremental:
                applycal(
                    vis=self.visfile,
                    spw=self.Imager.getSpw(),
                    spwmap=[self.spwmap, self.spwmap],
                    field=self.Imager.getField(),
                    gaintable=[self.input_caltable, caltable],
                    gainfield='',
                    calwt=False,
                    flagbackup=False,
                    interp=self.interp,
                    applymode=self.applymode
                )
            else:
                applycal(
                    vis=self.visfile,
                    spw=self.Imager.getSpw(),
                    spwmap=self.spwmap,
                    field=self.Imager.getField(),
                    gaintable=[caltable],
                    gainfield='',
                    calwt=False,
                    flagbackup=False,
                    interp=self.interp,
                    applymode=self.applymode
                )

            if self.flag_dataset_bool:
                self.flag_dataset(mode=self.flag_mode)

            imagename = self.imagename + '_ap' + str(i)

            self.Imager.run(imagename)

            self.psnr_history.append(self.Imager.getPSNR())

            print("Solint: " + str(self.solint[i]) + " - PSNR: " + str(self.psnr_history[i]))
            print("Noise: " + str(self.Imager.getSTDV() * 1000.0) + " mJy/beam")

            path_object = Path(self.visfile)
            current_visfile = "{0}_{2}{1}".format(
                Path.joinpath(path_object.parent, path_object.stem), path_object.suffix,
                "ap" + str(i)
            )
            if os.path.exists(current_visfile):
                shutil.rmtree(current_visfile)
            shutil.copytree(self.visfile, current_visfile)
            if self.restore_PSNR:
                if i > 0:
                    if self.psnr_history[-1] <= self.psnr_history[-2]:
                        self.restore_selfcal(caltable_version=self.caltables_versions[-1])
                        self.psnr_history.pop()
                        self.caltables.pop()
                        print(
                            "PSNR decreasing in this solution interval - restoring to last MS and exiting loop"
                        )
                        break
                    else:
                        print(
                            "PSNR improved on iteration {0} - Copying measurement set files...".
                            format(i)
                        )
                        if os.path.exists(current_visfile):
                            shutil.rmtree(current_visfile)
                        shutil.copytree(self.visfile, current_visfile)
                        self.visfile = current_visfile
                        self.Imager.inputvis = current_visfile
                        self.write_file_backup()

                elif self.selfcal_object:
                    if self.psnr_history[i] <= self.selfcal_object.getPSNRHistory()[-1]:
                        self.restore_selfcal(
                            caltable_version=self.selfcal_object.getCaltablesVersions()[-1]
                        )
                        self.psnr_history.pop()
                        self.caltables_versions.pop()
                        self.caltables.pop()
                        self.caltables = self.selfcal_object.getCaltables()
                        self.psnr_history = self.selfcal_object.getPSNRHistory()
                        self.caltables_versions = self.selfcal_object.getCaltablesVersions()
                        print(
                            "PSNR decreasing in this solution interval - restoring to last MS and exiting loop"
                        )
                        break
                    else:
                        print(
                            "PSNR improved on iteration {0} - Copying measurement set files...".
                            format(i)
                        )
                        if os.path.exists(current_visfile):
                            shutil.rmtree(current_visfile)
                        shutil.copytree(self.visfile, current_visfile)
                        self.visfile = current_visfile
                        self.Imager.inputvis = current_visfile
                        self.write_file_backup()
