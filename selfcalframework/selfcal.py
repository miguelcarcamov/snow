import os
import sys
import numpy as np
from flagmanager import flagmanager
from rmtables import rmtables
from gaincal import gaincal
from clearcal import clearcal
from delmod import delmod
from plotcal import plotcal
from applycal import applycal
from split import split
from flagdata import flagdata
import abc


class Selfcal(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, visfile="", Imager=None, refant="", spwmap=[], minblperant=4, want_plot=True, interp='linear', gaintype='T', solint=[], flag_mode="rflag", flag_dataset_bool=False, restore_PSNR=False, **kwargs):
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(**kwargs)
        self.caltables = []
        self.caltables_versions = []
        self.psnr_history = []

        if(self.Imager == None):
            print("Error, Imager Object is Nonetype")
            sys.exit(
                "Error, self-calibration objects cannot run without an imager object")

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

    def reset_selfcal(self, caltable_version=""):
        flagmanager(vis=self.visfile, mode='restore',
                    versionname=caltable_version)
        clearcal(self.visfile)
        delmod(self.visfile, otf=True)

    def restore_selfcal(self, caltable_version=""):
        flagmanager(vis=self.visfile, mode='restore',
                    versionname=caltable_version)
        delmod(self.visfile, otf=True)

    def flag_dataset(self, datacolumn="residuals", mode=""):
        flagdata(vis=self.visfile, datacolumn=datacolumn, action="apply", mode=mode, flagbackup=False)

    def plot_selfcal(self, caltable, xaxis="", yaxis="", timerange="", iteration="", antenna="", subplot=111, plotrange=[], want_plot=False, **kwargs):
        figfile_name = caltable + ".png"
        if want_plot:
            plotcal(caltable=caltable, xaxis=xaxis, yaxis=yaxis, timerange=timerange,
                    iteration=iteration, subplot=subplot, antenna=antenna, plotrange=plotrange, figfile=figfile_name, showgui=True)
        else:
            plotcal(caltable=caltable, xaxis=xaxis, yaxis=yaxis, timerange=timerange,
                    iteration=iteration, subplot=subplot, antenna=antenna, plotrange=plotrange, figfile=figfile_name, showgui=False)

    def selfcal_output(self, overwrite=False):
        outputvis = self.visfile + '.selfcal'

        if overwrite:
            os.system('rm -rf ' + outputvis)
        split(vis=self.visfile, outputvis=outputvis, datacolumn='corrected')
        return outputvis

    @abc.abstractmethod
    def run(self):
        return


class Ampcal(Selfcal):
    def __init__(self, visfile="", Imager=None, refant="", spwmap=[], minblperant=4, want_plot=True, interp='linear', gaintype='T', solint=[], flag_mode="rflag", flag_dataset_bool=False, restore_PSNR=False, selfcal_object=None, input_caltable="", **kwargs):

        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(**kwargs)

        super(Ampcal, self).__init__(visfile, Imager, refant, spwmap, minblperant, want_plot, interp, gaintype, solint, flag_mode, flag_dataset_bool, restore_PSNR, **kwargs)
        self.calmode = 'a'
        self.loops = len(self.solint)
        self.imagename = self.Imager.getOutputPath()

        if(self.selfcal_object == None and self.input_caltable == ""):
            print(
                "Error, Self-cal object is Nonetype and input_caltable is an empty string")
            sys.exit(
                "Error, Amplitude self-cal objects cannot run without an phase-cal object or input caltable")
        else:
            if(self.selfcal_object):
                self.input_caltable = self.selfcal_object.getCaltables()[-1]
            elif(self.input_caltable != ""):
                print("The caltable input must been already created")
                print("Self-cal table: " + self.input_caltable)
            else:
                print("Error, Ampcal needs a non-empty list of caltables")
                sys.exit(
                    "Error, Amplitude self-cal objects cannot run with an empty caltable list")

    def run(self):
        caltable = ""
        for i in range(0, self.loops):
            caltable = 'ampcal_' + str(i)
            self.caltables.append(caltable)
            rmtables(caltable)
            gaincal(vis=self.visfile, field=self.Imager.getField(), caltable=caltable, spw=self.Imager.getSpw(), gaintype=self.gaintype, refant=self.refant, calmode=self.calmode,
                    combine=self.combine, solint=self.solint[i], minsnr=self.minsnr, minblperant=self.minblperant, gaintable=self.input_caltable, spwmap=self.spwmap, solnorm=True)

            self.plot_selfcal(caltable, xaxis="time", yaxis="amp", iteration="antenna",
                              subplot=421, plotrange=[0, 0, 0.2, 1.8], want_plot=self.want_plot)

            versionname = 'before_ampcal_' + str(i)
            flagmanager(vis=self.visfile, mode='save',
                        versionname=versionname)
            self.caltables_versions.append(versionname)
            applycal(vis=self.visfile, spwmap=[self.spwmap, self.spwmap], field=self.Imager.getField(), gaintable=[
                     self.input_caltable, caltable], gainfield='', calwt=False, flagbackup=False, interp=self.interp)

            imagename = self.imagename + '_a' + str(i)

            self.Imager.run(imagename)

            if(self.flag_dataset_bool):
                flag_dataset(mode=self.flag_mode)
            self.psnr_history.append(self.Imager.getPSNR())

            print("Solint: " + str(self.solint[i]) +
                  " - PSNR: " + str(self.psnr_history[i]))
            print("Noise: " + str(self.Imager.getSTDV() * 1000.0) + " mJy/beam")
            if(self.restore_PSNR):
                if(i > 0):
                    if(self.psnr_history[i] < self.psnr_history[i - 1]):
                        self.restore_selfcal(
                            caltable_version=self.caltables_versions[i - 1])
                        self.psnr_history.pop()
                        self.caltables_versions.pop()
                        self.caltables.pop()
                        print(
                            "PSNR decreasing in this solution interval - restoring to last MS and exiting loop")
                        break
                else:
                    if(self.selfcal_object):
                        if(self.psnr_history[i] < self.selfcal_object.getPSNRHistory()[-1]):
                            self.restore_selfcal(
                                caltable_version=self.selfcal_object.getCaltablesVersions()[-1])
                            self.psnr_history.pop()
                            self.caltables_versions.pop()
                            self.caltables.pop()
                            print(
                                "PSNR decreasing in this solution interval - restoring to last MS and exiting loop")
                            break


class Phasecal(Selfcal):
    def __init__(self, visfile="", Imager=None, refant="", spwmap=[], minblperant=4, want_plot=True, interp='linear', gaintype='T', solint=[], flag_mode="rflag", flag_dataset_bool=False, restore_PSNR=False, input_caltable="", **kwargs):

        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(**kwargs)
        super(Phasecal, self).__init__(visfile, Imager, refant, spwmap, minblperant, want_plot, interp, gaintype, solint, flag_mode, flag_dataset_bool, restore_PSNR, **kwargs)

        self.calmode = 'p'
        self.loops = len(self.solint)
        self.imagename = self.Imager.getOutputPath()

    def run(self):
        caltable = "before_selfcal"
        flagmanager(vis=self.visfile, mode='save',
                    versionname=caltable)
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

            gaincal(vis=self.visfile, caltable=caltable, field=self.Imager.getField(), spw=self.Imager.getSpw(), gaintype=self.gaintype, refant=self.refant,
                    calmode=self.calmode, combine=self.combine, solint=self.solint[i], minsnr=self.minsnr, minblperant=self.minblperant)

            self.plot_selfcal(caltable, xaxis="time", yaxis="phase", iteration="antenna",
                              subplot=421, plotrange=[0, 0, -180, 180], want_plot=self.want_plot)

            versionname = 'before_phasecal_' + str(i)
            flagmanager(vis=self.visfile, mode='save',
                        versionname=versionname)
            self.caltables_versions.append(versionname)

            applycal(vis=self.visfile, field=self.Imager.getField(), spwmap=self.spwmap, gaintable=[
                     caltable], gainfield='', calwt=False, flagbackup=False, interp=self.interp)

            imagename = self.imagename + '_ph' + str(i)

            self.Imager.run(imagename)

            if(self.flag_dataset_bool):
                flag_dataset(mode=self.flag_mode)

            self.psnr_history.append(self.Imager.getPSNR())
            print("Solint: " + str(self.solint[i]) +
                  " - PSNR: " + str(self.psnr_history[-1]))
            print("Noise: " + str(self.Imager.getSTDV() * 1000.0) + " mJy/beam")
            if(self.restore_PSNR):
                if(self.psnr_history[-1] < self.psnr_history[-2]):
                    self.restore_selfcal(
                        caltable_version=self.caltables_versions[i])
                    self.psnr_history.pop()
                    self.caltables_versions.pop()
                    self.caltables.pop()
                    print(
                        "PSNR decreasing in this solution interval - restoring to last MS and exiting loop")
                    break


class AmpPhasecal(Selfcal):
    def __init__(self, visfile="", Imager=None, refant="", spwmap=[], minblperant=4, want_plot=True, interp='linear', gaintype='T', solint=[], flag_mode="rflag", flag_dataset_bool=False, restore_PSNR=False, selfcal_object=None, **kwargs):

        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(**kwargs)
        super(AmpPhasecal, self).__init__(visfile, Imager, refant, spwmap, minblperant, want_plot, interp, gaintype, solint, flag_mode, flag_dataset_bool, restore_PSNR, **kwargs)

        self.calmode = 'ap'
        self.loops = len(self.solint)
        self.imagename = self.Imager.getOutputPath()

        if(self.selfcal_object == None and self.input_caltable == ""):
            print(
                "Error, Self-cal object is Nonetype and input_caltable is an empty string")
            sys.exit(
                "Error, Amplitude self-cal objects cannot run without an phase-cal object or input caltable")
        else:
            if(self.selfcal_object):
                self.input_caltable = self.selfcal_object.getCaltables()[-1]
            elif(self.input_caltable != ""):
                print("The caltable input must been already created")
                print("Self-cal table: " + self.input_caltable)
            else:
                print("Error, Ampcal needs a non-empty list of caltables")
                sys.exit(
                    "Error, Amplitude self-cal objects cannot run with an empty caltable list")

    def run(self):
        caltable = ""
        for i in range(0, self.loops):
            caltable = 'apcal_' + str(i)
            self.caltables.append(caltable)
            rmtables(caltable)
            gaincal(vis=self.visfile, field=self.Imager.getField(), caltable=caltable, spw=self.Imager.getSpw(), gaintype=self.gaintype, refant=self.refant, calmode=self.calmode,
                    combine=self.combine, solint=self.solint[
                        i], minsnr=self.minsnr, minblperant=self.minblperant, gaintable=self.input_caltable, spwmap=self.spwmap,
                    solnorm=True)

            self.plot_selfcal(caltable, xaxis="time", yaxis="amp", iteration="antenna",
                              subplot=421, plotrange=[0, 0, 0.2, 1.8], want_plot=self.want_plot)

            versionname = 'before_apcal_' + str(i)
            flagmanager(vis=self.visfile, mode='save',
                        versionname=versionname)
            self.caltables_versions.append(versionname)

            applycal(vis=self.visfile, spwmap=[self.spwmap, self.spwmap], field=self.Imager.getField(), gaintable=[
                     self.input_caltable, caltable], gainfield='', calwt=False, flagbackup=False, interp=self.interp)

            imagename = self.imagename + '_ap' + str(i)

            self.Imager.run(imagename)

            if(self.flag_dataset_bool):
                flag_dataset(mode=self.flag_mode)

            self.psnr_history.append(self.Imager.getPSNR())
            print("Solint: " + str(self.solint[i]) +
                  " - PSNR: " + str(self.psnr_history[i]))
            print("Noise: " + str(self.Imager.getSTDV() * 1000.0) + " mJy/beam")
            if(self.restore_PSNR):
                if(i > 0):
                    if(self.psnr_history[i] < self.psnr_history[i - 1]):
                        self.restore_selfcal(
                            caltable_version=self.caltables_versions[i - 1])
                        self.psnr_history.pop()
                        self.caltables_versions.pop()
                        self.caltables.pop()
                        print(
                            "PSNR decreasing in this solution interval - restoring to last MS and exiting loop")
                        break
                else:
                    if(self.selfcal_object):
                        if(self.psnr_history[i] < self.selfcal_object.getPSNRHistory()[-1]):
                            self.restore_selfcal(
                                caltable_version=self.selfcal_object.getCaltablesVersions()[-1])
                            self.psnr_history.pop()
                            self.caltables_versions.pop()
                            self.caltables.pop()
                            print(
                                "PSNR decreasing in this solution interval - restoring to last MS and exiting loop")
                            break
