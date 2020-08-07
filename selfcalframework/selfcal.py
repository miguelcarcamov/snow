import os
import numpy as np
from flagmanager import flagmanager
from rmtables import rmtables
from gaincal import gaincal
from plotcal import plotcal
from applycal import applycal
from split import split


class Selfcal(object):

    def __init__(self, visfile="", Imager=None, refant="", spwmap=[], minblperant=4, want_plot=True, interp='linear', **kwargs):
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(kwargs)

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
    def __init__(self, solint=[], selfcal_object=None, **kwargs):
        super(Ampcal, self).__init__(selfcal_object.visfile, selfcal_object.Imager, selfcal_object.refant,
                                     selfcal_object.spwmap, selfcal_object.minblperant, selfcal_object.want_plot, selfcal_object.interp, **kwargs)

        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(kwargs)
        self.calmode = 'a'
        self.loops = len(self.solint)
        self.visfile = selfcal_object.visfile
        self.Imager = selfcal_object.Imager
        self.minblperant = selfcal_object.minblperant
        self.refant = selfcal_object.refant
        self.spwmap = selfcal_object.spwmap
        self.interp = selfcal_object.interp
        self.want_plot = selfcal_object.want_plot
        self.imagename = self.Imager.getOutputPath()

    def run(self):
        caltable = ""
        for i in range(0, self.loops):
            caltable = 'ampcal_' + str(i)
            rmtables(caltable)
            gaincal(vis=self.visfile, field=self.Imager.getField(), caltable=caltable, spw=self.Imager.getSpw(), gaintype='G', refant=self.refant, calmode=self.calmode,
                    combine=self.combine, solint=self.solint[i], minsnr=self.minsnr, minblperant=self.minblperant, gaintable=self.input_caltable, spwmap=self.spwmap, solnorm=True)

            self.plot_selfcal(caltable, xaxis="time", yaxis="amp", iteration="antenna",
                              subplot=421, plotrange=[0, 0, 0.2, 1.8], want_plot=self.want_plot)

            applycal(vis=self.visfile, spwmap=[self.spwmap, self.spwmap], field=self.Imager.getField(), gaintable=[
                     self.input_caltable, caltable], gainfield='', calwt=False, flagbackup=False, interp=self.interp)

            flagmanager(vis=self.visfile, mode='save',
                        versionname='after_ampcal' + str(i))

            imagename = self.imagename + '_a' + str(i)

            self.Imager.run(imagename)
        return caltable


class Phasecal(Selfcal):
    def __init__(self, solint=[], selfcal_object=None, **kwargs):
        super(Phasecal, self).__init__(selfcal_object.visfile, selfcal_object.Imager, selfcal_object.refant,
                                       selfcal_object.spwmap, selfcal_object.minblperant, selfcal_object.want_plot, selfcal_object.interp, **kwargs)
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(kwargs)
        self.calmode = 'p'
        self.loops = len(self.solint)
        self.visfile = selfcal_object.visfile
        self.Imager = selfcal_object.Imager
        self.minblperant = selfcal_object.minblperant
        self.refant = selfcal_object.refant
        self.spwmap = selfcal_object.spwmap
        self.want_plot = selfcal_object.want_plot
        self.interp = selfcal_object.interp
        self.imagename = self.Imager.getOutputPath()

    def run(self):
        flagmanager(vis=self.visfile, mode='save',
                    versionname='before_phasecal', merge='replace')
        caltable = ""
        for i in range(0, self.loops):
            imagename = self.imagename + '_ph' + str(i)
            self.Imager.run(imagename)
            caltable = 'pcal' + str(i)
            rmtables(caltable)

            gaincal(vis=self.visfile, caltable=caltable, field=self.Imager.getField(), spw=self.Imager.getSpw(), gaintype='G', refant=self.refant,
                    calmode=self.calmode, combine=self.combine, solint=self.solint[i], minsnr=self.minsnr, minblperant=self.minblperant)

            self.plot_selfcal(caltable, xaxis="time", yaxis="phase", iteration="antenna",
                              subplot=421, plotrange=[0, 0, -180, 180], want_plot=self.want_plot)

            applycal(vis=self.visfile, field=self.Imager.getField(), spwmap=self.spwmap, gaintable=[
                     caltable], gainfield='', calwt=False, flagbackup=False, interp=self.interp)

            flagmanager(vis=self.visfile, mode='save',
                        versionname='after_' + caltable)
        return caltable


class AmpPhasecal(Selfcal):
    def __init__(self, solint=[], selfcal_object=None, **kwargs):
        super(AmpPhasecal, self).__init__(selfcal_object.visfile, selfcal_object.Imager, selfcal_object.refant,
                                          selfcal_object.spwmap, selfcal_object.minblperant, selfcal_object.want_plot, selfcal_object.interp, **kwargs)

        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(kwargs)
        self.calmode = 'ap'
        self.loops = len(self.solint)
        self.visfile = selfcal_object.visfile
        self.Imager = selfcal_object.Imager
        self.minblperant = selfcal_object.minblperant
        self.refant = selfcal_object.refant
        self.spwmap = selfcal_object.spwmap
        self.want_plot = selfcal_object.want_plot
        self.interp = selfcal_object.interp
        self.imagename = self.Imager.getOutputPath()

    def run(self):
        caltable = ""
        for i in range(0, self.loops):
            caltable = 'apcal_' + str(i)
            rmtables(caltable)
            gaincal(vis=self.visfile, field=self.Imager.getField(), caltable=caltable, spw=self.Imager.getSpw(), gaintype='G', refant=self.refant, calmode=self.calmode,
                    combine=self.combine, solint=self.solint[i], minsnr=self.minsnr, minblperant=self.minblperant, gaintable=self.input_caltable, spwmap=self.spwmap,
                    solnorm=True)

            self.plot_selfcal(caltable, xaxis="time", yaxis="amp", iteration="antenna",
                              subplot=421, plotrange=[0, 0, 0.2, 1.8], want_plot=self.want_plot)

            applycal(vis=self.visfile, spwmap=[self.spwmap, self.spwmap], field=self.Imager.getField(), gaintable=[
                     self.input_caltable, caltable], gainfield='', calwt=False, flagbackup=False, interp=self.interp)

            flagmanager(vis=self.visfile, mode='save',
                        versionname='after_apcal' + str(i))

            imagename = self.imagename + '_ap' + str(i)

            self.Imager.run(imagename)
