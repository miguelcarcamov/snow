import numpy as np
import os
import time
from tclean import tclean
from image_utils import *
import abc


class Imager(object):
    __metaclass__ = abc.ABCMeta
    def __init__(self, robust=2.0, field="", spw="", savemodel=True, verbose=True, **kwargs):
        self.psnr = 0.0
        self.peak = 0.0
        self.stdv = 0.0
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(**kwargs)

    def getField(self):
        return self.field

    def getSpw(self):
        return self.spw

    def getVis(self):
        return self.inputvis

    def getOutputPath(self):
        return self.output

    def getPSNR(self):
        return self.psnr

    def getPeak(self):
        return self.peak

    def getSTDV(self):
        return self.stdv

    def calculateStatistics_FITS(self, signal_fits_name="", residual_fits_name="", stdv_pixels=80):
        self.psnr, self.peak, self.stdv = calculatePSNR_FITS(
            signal_fits_name, residual_fits_name, stdv_pixels)

    def calculateStatistics_MSImage(self, signal_ms_name="", residual_ms_name="", stdv_pixels=80):
        self.psnr, self.peak, self.stdv = calculatePSNR_MS(
            signal_ms_name, residual_ms_name, stdv_pixels)

    @abc.abstractmethod
    def run(self, imagename=""):
        return


class Clean(Imager):
    def __init__(self, nterms=1, threshold=0.0, nsigma=0.0, interactive=False, usemask="auto-multithresh", negativethreshold=0.0, lownoisethreshold=1.5, noisethreshold=4.25,
                 sidelobethreshold=2.0, minbeamfrac=0.3, deconvolver="hogbom", scales=[], pbcor=False, cycleniter=-1, savemodel=True, clean_savemodel=None, **kwargs):
        super(Clean, self).__init__(**kwargs)
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(**kwargs)

        if(self.savemodel):
            self.clean_savemodel = "modelcolumn"

    def run(self, imagename=""):
        imsize = [self.M, self.N]
        tclean(vis=self.inputvis, imagename=imagename, field=self.field,
               datacolumn=self.datacolumn, specmode=self.specmode, stokes=self.stokes, deconvolver=self.deconvolver, scales=self.scales, nterms=self.nterms,
               imsize=imsize, cell=self.cell, weighting="briggs", robust=self.robust, niter=self.niter, threshold=self.threshold, nsigma=self.nsigma,
               interactive=self.interactive, gridder=self.gridder, pbcor=self.pbcor, savemodel=self.clean_savemodel, usemask=self.usemask,
               negativethreshold=self.negativethreshold, lownoisethreshold=self.lownoisethreshold, noisethreshold=self.noisethreshold,
               sidelobethreshold=self.sidelobethreshold, minbeamfrac=self.minbeamfrac, cycleniter=self.cycleniter, verbose=self.verbose)

        if(self.deconvolver != "mtmfs"):
            restored_image = imagename + ".image"
            residual_image = imagename + ".residual"
        else:
            restored_image = imagename + ".image.tt0"
            residual_image = imagename + ".residual.tt0"

        self.calculateStatistics_MSImage(
            signal_ms_name=restored_image, residual_ms_name=residual_image)


class WSClean(Imager):
    def __init__(self, **kwargs):
        super(Clean, self).__init__(**kwargs)
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(**kwargs)


class GPUvmem(Imager):
    def __init(self, **kwargs):
        super(GPUvmem, self).__unut__(**kwargs)
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(**kwargs)
