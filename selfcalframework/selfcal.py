import os
import sys
import shutil
import numpy as np
from casatasks import flagmanager
from casatasks import rmtables
from casatasks import gaincal
from casatasks import clearcal
from casatasks import delmod
from casatasks import applycal
from casatasks import split
from casatasks import flagdata
from casatasks import uvsub
from casatools import table
from casatasks import statwt
from abc import ABCMeta, abstractmethod

tb = table()


class Selfcal(metaclass=ABCMeta):

    def __init__(self, visfile="", Imager=None, refant="", spwmap=[], minblperant=4, want_plot=True, interp='linear',
                 gaintype='T', uvrange="", solint=[], varchange=None, minsnr=3.0, applymode="calflag",
                 flag_mode="rflag", combine="", flag_dataset_bool=False, restore_PSNR=False, subtract_source=False):
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.caltables = []
        self.caltables_versions = []
        self.psnr_history = []

        if self.Imager is None:
            print("Error, Imager Object is Nonetype")
            sys.exit(
                "Error, self-calibration objects cannot run without an imager object")

        if self.varchange is not None:
            list_of_values = [value for key, value in self.varchange.items()]
            it = iter(list_of_values)
            if not all(len(l) == len(self.solint) for l in it):
                sys.exit("Error, length of solint and variable that changes through iterations must be the same")

        if self.subtract_source == True:
            if self.Imager.getPhaseCenter() != "":
                print("Error, phase center needs to be set if a source is going to be subtracted")

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

    def save_selfcal(self, caltable_version="", overwrite=True):
        if overwrite:
            flagmanager(vis=self.visfile, mode='delete', versionname=caltable_version)
        flagmanager(vis=self.visfile, mode='save', versionname=caltable_version)

    def reset_selfcal(self, caltable_version=""):
        flagmanager(vis=self.visfile, mode='restore',
                    versionname=caltable_version)
        clearcal(self.visfile)
        delmod(vis=self.visfile, otf=True, scr=True)

    def restore_selfcal(self, caltable_version=""):
        flagmanager(vis=self.visfile, mode='restore',
                    versionname=caltable_version)
        delmod(vis=self.visfile, otf=True, scr=True)

    def flag_dataset(self, datacolumn="residual", mode="rflag"):
        flagdata(vis=self.visfile, mode=mode, datacolumn=datacolumn, field='', timecutoff=5.0,
                 freqcutoff=5.0, timefit='line', freqfit='line', flagdimension='freqtime',
                 extendflags=False, timedevscale=3.0, freqdevscale=3.0, spectralmax=500,
                 extendpols=False, growaround=False, flagneartime=False, flagnearfreq=False,
                 action='apply', flagbackup=True, overwrite=True, writeflags=True)

    def ismodel_in_dataset(self):
        tb.open(tablename=self.visfile)
        columns = tb.colnames()
        tb.close()
        if "MODEL_DATA" in columns:
            return True
        else:
            return False

    def set_Imager_attributes_from_dict(self, iteration=0):
        if self.varchange is not None:
            for key in self.varchange.keys():
                setattr(self.Imager, key, self.varchange[key][iteration])

    def plot_selfcal(self, caltable, xaxis="", yaxis="", iteration="", timerange="", antenna="", subplot=[1, 1],
                     plotrange=[], want_plot=False, **kwargs):
        figfile_name = caltable + ".png"
        if want_plot:
            print(
                "Plot selfcal is trying to plot, but the function is not part of CASA 6 yet")
            # au.plotms(vis=caltable, xaxis=xaxis, yaxis=yaxis, iteration=iteration, gridrows=subplot[0],
            # gridcols=subplot[1], antenna=antenna, timerange=timerange, plotrange=plotrange, plotfile=figfile_name,
            # overwrite=True, showgui=True)
        else:
            print(
                "Plot selfcal is trying to plot, but the function is not part of CASA 6 yet")
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

        if self.selfcal_object is None and self.input_caltable == "":
            print(
                "Error, Self-cal object is Nonetype and input_caltable is an empty string")
            sys.exit(
                "Error, Amplitude self-cal objects cannot run without an phase-cal object or input caltable")
        else:
            if self.selfcal_object:
                self.input_caltable = self.selfcal_object.getCaltables()[-1]
            elif self.input_caltable != "":
                if not os.path.exists(self.input_caltable):
                    sys.exit("The caltable " + self.input_caltable + " needs to be created")
            else:
                print("Error, Ampcal needs a non-empty list of caltables")
                sys.exit(
                    "Error, Amplitude self-cal objects cannot run with an empty caltable list")

    def run(self):
        caltable = ""
        if not self.ismodel_in_dataset():
            imagename = self.imagename + "before_apcal"
            self.Imager.run(imagename)
            print("Before amplitude self-cal: - PSNR: " + str(self.Imager.getPSNR()))
            print("Noise: " + str(self.Imager.getSTDV() * 1000.0) + " mJy/beam")
            self.psnr_history.append(self.Imager.getPSNR())

        for i in range(0, self.loops):
            caltable = 'ampcal_' + str(i)
            self.caltables.append(caltable)
            rmtables(caltable)
            gaincal(vis=self.visfile, field=self.Imager.getField(), caltable=caltable, spw=self.Imager.getSpw(),
                    uvrange=self.uvrange, gaintype=self.gaintype, refant=self.refant, calmode=self.calmode,
                    combine=self.combine, solint=self.solint[i], minsnr=self.minsnr, minblperant=self.minblperant,
                    gaintable=self.input_caltable, spwmap=self.spwmap, solnorm=True)

            self.plot_selfcal(caltable, xaxis="time", yaxis="amp", iteration="antenna",
                              subplot=[4, 2], plotrange=[0, 0, 0.2, 1.8], want_plot=self.want_plot)

            versionname = 'before_ampcal_' + str(i)
            self.save_selfcal(caltable_version=versionname, overwrite=True)
            self.caltables_versions.append(versionname)
            applycal(vis=self.visfile, spw=self.Imager.getSpw(), spwmap=[self.spwmap, self.spwmap],
                     field=self.Imager.getField(), gaintable=[
                    self.input_caltable, caltable], gainfield='', calwt=False, flagbackup=False, interp=self.interp,
                     applymode=self.applymode)

            if self.flag_dataset_bool:
                self.flag_dataset(mode=self.flag_mode)

            imagename = self.imagename + '_a' + str(i)

            self.set_Imager_attributes_from_dict(i)

            self.Imager.run(imagename)

            self.psnr_history.append(self.Imager.getPSNR())

            print("Solint: " + str(self.solint[i]) +
                  " - PSNR: " + str(self.psnr_history[i]))
            print("Noise: " + str(self.Imager.getSTDV() * 1000.0) + " mJy/beam")
            if self.restore_PSNR:
                if i > 0:
                    if self.psnr_history[i] <= self.psnr_history[i - 1]:
                        self.restore_selfcal(
                            caltable_version=self.caltables_versions[i - 1])
                        self.psnr_history.pop()
                        self.caltables_versions.pop()
                        self.caltables.pop()
                        print(
                            "PSNR decreasing in this solution interval - restoring to last MS and exiting loop")
                        break
                else:
                    if self.selfcal_object:
                        if self.psnr_history[i] <= self.selfcal_object.getPSNRHistory()[-1]:
                            self.restore_selfcal(
                                caltable_version=self.selfcal_object.getCaltablesVersions()[-1])
                            self.psnr_history.pop()
                            self.caltables_versions.pop()
                            self.caltables.pop()
                            self.caltables = self.selfcal_object.getCaltables()
                            self.psnr_history = self.selfcal_object.getPSNRHistory()
                            self.caltables_versions = self.selfcal_object.getCaltablesVersions()
                            print(
                                "PSNR decreasing or equal in this solution interval - restoring to last MS and "
                                "exiting loop")
                            break


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

    def run(self):
        caltable = "before_selfcal"
        self.save_selfcal(caltable_version=caltable, overwrite=True)
        self.caltables_versions.append(caltable)
        imagename = self.imagename + '_original'
        self.Imager.run(imagename)
        print("Original: - PSNR: " + str(self.Imager.getPSNR()))
        print("Noise: " + str(self.Imager.getSTDV() * 1000.0) + " mJy/beam")
        self.psnr_history.append(self.Imager.getPSNR())

        for i in range(0, self.loops):
            caltable = 'pcal' + str(i)
            self.caltables.append(caltable)
            rmtables(caltable)

            gaincal(vis=self.visfile, caltable=caltable, field=self.Imager.getField(), spw=self.Imager.getSpw(),
                    uvrange=self.uvrange, gaintype=self.gaintype, refant=self.refant,
                    calmode=self.calmode, combine=self.combine, solint=self.solint[i], minsnr=self.minsnr,
                    spwmap=self.spwmap, minblperant=self.minblperant)

            self.plot_selfcal(caltable, xaxis="time", yaxis="phase", iteration="antenna",
                              subplot=[4, 2], plotrange=[0, 0, -180, 180], want_plot=self.want_plot)

            versionname = 'before_phasecal_' + str(i)
            self.save_selfcal(caltable_version=versionname, overwrite=True)
            self.caltables_versions.append(versionname)

            applycal(vis=self.visfile, field=self.Imager.getField(), spw=self.Imager.getSpw(), spwmap=self.spwmap,
                     gaintable=[
                         caltable], gainfield='', calwt=False, flagbackup=False, interp=self.interp,
                     applymode=self.applymode)

            if self.flag_dataset_bool:
                self.flag_dataset(mode=self.flag_mode)

            imagename = self.imagename + '_ph' + str(i)

            self.set_Imager_attributes_from_dict(i)

            self.Imager.run(imagename)

            self.psnr_history.append(self.Imager.getPSNR())

            print("Solint: " + str(self.solint[i]) +
                  " - PSNR: " + str(self.psnr_history[-1]))
            print("Noise: " + str(self.Imager.getSTDV() * 1000.0) + " mJy/beam")
            if self.restore_PSNR:
                if self.psnr_history[-1] <= self.psnr_history[-2]:
                    self.restore_selfcal(
                        caltable_version=self.caltables_versions[i])
                    self.psnr_history.pop()
                    self.caltables_versions.pop()
                    self.caltables.pop()
                    print(
                        "PSNR decreasing or equal in this solution interval - restoring to last MS and exiting loop")
                    break


class AmpPhasecal(Selfcal):
    def __init__(self, selfcal_object=None, input_caltable="", **kwargs):
        super(AmpPhasecal, self).__init__(**kwargs)
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])

        self.calmode = 'ap'
        self.loops = len(self.solint)
        self.imagename = self.Imager.getOutputPath()

        if self.selfcal_object == None and self.input_caltable == "":
            print(
                "Error, Self-cal object is Nonetype and input_caltable is an empty string")
            sys.exit(
                "Error, Amplitude self-cal objects cannot run without an phase-cal object or input caltable")
        else:
            if self.selfcal_object:
                self.input_caltable = self.selfcal_object.getCaltables()[-1]
            elif self.input_caltable != "":
                if not os.path.exists(self.input_caltable):
                    sys.exit("The caltable " + self.input_caltable + " needs to be created")
            else:
                print("Error, Ampcal needs a non-empty list of caltables")
                sys.exit(
                    "Error, Amplitude self-cal objects cannot run with an empty caltable list")

    def run(self):
        caltable = ""
        if not self.ismodel_in_dataset():
            imagename = self.imagename + "before_apcal"
            self.Imager.run(imagename)
            print("Before amplitude-phase self-cal: - PSNR: " + str(self.Imager.getPSNR()))
            print("Noise: " + str(self.Imager.getSTDV() * 1000.0) + " mJy/beam")
            self.psnr_history.append(self.Imager.getPSNR())

        for i in range(0, self.loops):
            caltable = 'apcal_' + str(i)
            self.caltables.append(caltable)
            rmtables(caltable)
            gaincal(vis=self.visfile, field=self.Imager.getField(), caltable=caltable, spw=self.Imager.getSpw(),
                    uvrange=self.uvrange, gaintype=self.gaintype, refant=self.refant, calmode=self.calmode,
                    combine=self.combine, solint=self.solint[
                    i], minsnr=self.minsnr, minblperant=self.minblperant, gaintable=self.input_caltable,
                    spwmap=self.spwmap,
                    solnorm=True)

            self.plot_selfcal(caltable, xaxis="time", yaxis="amp", iteration="antenna",
                              subplot=[4, 2], plotrange=[0, 0, 0.2, 1.8], want_plot=self.want_plot)

            versionname = 'before_apcal_' + str(i)
            self.save_selfcal(caltable_version=versionname, overwrite=True)
            self.caltables_versions.append(versionname)

            applycal(vis=self.visfile, spw=self.Imager.getSpw(), spwmap=[self.spwmap, self.spwmap],
                     field=self.Imager.getField(), gaintable=[
                    self.input_caltable, caltable], gainfield='', calwt=False, flagbackup=False, interp=self.interp,
                     applymode=self.applymode)

            if self.flag_dataset_bool:
                self.flag_dataset(mode=self.flag_mode)

            imagename = self.imagename + '_ap' + str(i)

            self.set_Imager_attributes_from_dict(i)

            self.Imager.run(imagename)

            self.psnr_history.append(self.Imager.getPSNR())

            print("Solint: " + str(self.solint[i]) +
                  " - PSNR: " + str(self.psnr_history[i]))
            print("Noise: " + str(self.Imager.getSTDV() * 1000.0) + " mJy/beam")
            if self.restore_PSNR:
                if i > 0:
                    if self.psnr_history[i] <= self.psnr_history[i - 1]:
                        self.restore_selfcal(
                            caltable_version=self.caltables_versions[i - 1])
                        self.psnr_history.pop()
                        self.caltables_versions.pop()
                        self.caltables.pop()
                        print(
                            "PSNR decreasing in this solution interval - restoring to last MS and exiting loop")
                        break
                else:
                    if self.selfcal_object:
                        if self.psnr_history[i] <= self.selfcal_object.getPSNRHistory()[-1]:
                            self.restore_selfcal(
                                caltable_version=self.selfcal_object.getCaltablesVersions()[-1])
                            self.psnr_history.pop()
                            self.caltables_versions.pop()
                            self.caltables.pop()
                            self.caltables = self.selfcal_object.getCaltables()
                            self.psnr_history = self.selfcal_object.getPSNRHistory()
                            self.caltables_versions = self.selfcal_object.getCaltablesVersions()
                            print(
                                "PSNR decreasing in this solution interval - restoring to last MS and exiting loop")
                            break
