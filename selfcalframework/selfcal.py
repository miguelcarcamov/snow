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


class Selfcal(object):

    def __init__(self, visfile="", Imager=None, refant="", spwmap=[], minblperant=4, want_plot=True, interp='linear', gaintype='T', restore_PSNR=False, **kwargs):
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
        return psnr_history

    def reset_selfcal(self, caltable_version=""):
        flagmanager(vis=self.visfile, mode='restore',
                    versionname=caltable_version)
        clearcal(self.visfile)
        delmod(self.visfile, otf=True)

    def restore_selfcal(self, caltable_version=""):
        flagmanager(vis=self.visfile, mode='restore',
                    versionname=caltable_version)
        delmod(self.visfile, otf=True)

    def plot_selfcal(self, caltable, xaxis="", yaxis="", timerange="", iteration="", antenna="", subplot=111, plotrange=[], want_plot=False, **kwargs):
        if want_plot:
            plotcal(caltable=caltable, xaxis=xaxis, yaxis=yaxis, timerange=timerange,
                    iteration=iteration, subplot=subplot, antenna=antenna, plotrange=plotrange)

    def selfcal_output(self, overwrite=False):
        outputvis = self.visfile + '.selfcal'

        if overwrite:
            os.system('rm -rf ' + outputvis)
        split(vis=self.visfile, outputvis=outputvis, datacolumn='corrected')
        return outputvis


class Ampcal(Selfcal):
    def __init__(self, visfile="", Imager=None, selfcal_object=None, refant="", spwmap=[], minblperant=4, want_plot=True, interp='linear', gaintype='T', solint=[], **kwargs):

        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(**kwargs)

        super(Ampcal, self).__init__(visfile, Imager, refant, spwmap,
                                     minblperant, want_plot, interp, gaintype, **kwargs)
        self.calmode = 'a'
        self.loops = len(self.solint)
        self.imagename = self.Imager.getOutputPath()

        if(self.selfcal_object == None):
            print("Error, Self-cal object is Nonetype")
            sys.exit(
                "Error, Amplitude self-cal objects cannot run without an phase-cal object")

    def run(self):
        caltable = ""
        input_caltable = self.selfcal_object.getCaltables()[-1]
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
                     input_caltable, caltable], gainfield='', calwt=False, flagbackup=False, interp=self.interp)

            imagename = self.imagename + '_a' + str(i)

            self.Imager.run(imagename)

            self.psnr_history.append(self.Imager.getPSNR())

            print("Solint: " + self.solint[i] +
                  " - PSNR: " + self.psnr_history[i])
            if(self.restore_PSNR):
                if(i > 0):
                    if(self.psnr_history[i] < self.psnr_history[i - 1]):
                        self.restore_selfcal(
                            caltable_version=self.caltables_versions[i - 1])
                        sys.exit(
                            "PSNR decreasing in this solution interval - restoring to last MS")
                else:
                    if(self.psnr_history[i] < self.selfcal_object.getPSNRHistory()[-1]):
                        self.restore_selfcal(
                            caltable_version=self.selfcal_object.getCaltablesVersions()[-1])
                        sys.exit(
                            "PSNR decreasing in this solution interval - restoring to last MS")


class Phasecal(Selfcal):
    def __init__(self, visfile="", Imager=None, refant="", spwmap=[], minblperant=4, want_plot=True, interp='linear', gaintype='T', solint=[], **kwargs):

        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(**kwargs)
        super(Phasecal, self).__init__(visfile, Imager, refant, spwmap,
                                       minblperant, want_plot, interp, gaintype, **kwargs)

        self.calmode = 'p'
        self.loops = len(self.solint)
        self.imagename = self.Imager.getOutputPath()

    def run(self):
        caltable = ""
        flagmanager(vis=self.visfile, mode='save',
                    versionname="before_selfcal")
        for i in range(0, self.loops):
            imagename = self.imagename + '_ph' + str(i)
            self.Imager.run(imagename)
            self.psnr_history.append(self.Imager.getPSNR())
            print("Solint: " + self.solint[i] +
                  " - PSNR: " + self.psnr_history[i])

            if(self.restore_PSNR):
                if(i > 0):
                    if(self.psnr_history[i] < self.psnr_history[i - 1]):
                        self.restore_selfcal(
                            caltable_version=self.caltables_versions[i - 1])
                        sys.exit(
                            "PSNR decreasing in this solution interval - restoring to last MS")

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


class AmpPhasecal(Selfcal):
    def __init__(self, visfile="", Imager=None, selfcal_object=None, refant="", spwmap=[], minblperant=4, want_plot=True, interp='linear', gaintype='T', solint=[], **kwargs):

        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(**kwargs)
        super(AmpPhasecal, self).__init__(visfile, Imager, refant,
                                          spwmap, minblperant, want_plot, interp, gaintype, **kwargs)

        self.calmode = 'ap'
        self.loops = len(self.solint)
        self.imagename = self.Imager.getOutputPath()

        if(self.selfcal_object == None):
            print("Error, Self-cal object is Nonetype")
            sys.exit(
                "Error, Amplitude-phase self-cal objects cannot run without a phase-cal/amplitude-cal object")

    def run(self):
        caltable = ""
        input_caltable = self.selfcal_object.getCaltables()[-1]
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
                     input_caltable, caltable], gainfield='', calwt=False, flagbackup=False, interp=self.interp)

            imagename = self.imagename + '_ap' + str(i)

            self.Imager.run(imagename)
            self.psnr_history.append(self.Imager.getPSNR())
            print("Solint: " + self.solint[i] +
                  " - PSNR: " + self.psnr_history[i])
            if(self.restore_PSNR):
                if(i > 0):
                    if(self.psnr_history[i] < self.psnr_history[i - 1]):
                        self.restore_selfcal(
                            caltable_version=self.caltables_versions[i - 1])
                        sys.exit(
                            "PSNR decreasing in this solution interval - restoring to last MS")
                else:
                    if(self.psnr_history[i] < self.selfcal_object.getPSNRHistory()[-1]):
                        self.restore_selfcal(
                            caltable_version=self.selfcal_object.getCaltablesVersions()[-1])
                        sys.exit(
                            "PSNR decreasing in this solution interval - restoring to last MS")
