import numpy as np
import os
import time
from tclean import tclean
from importfits import importfits
from exportfits import exportfits
from imhead import imhead
from immath import immath
from qa import qa
from ia import ia
from image_utils import *
import abc
import subprocess


class Imager(object):
    __metaclass__ = abc.ABCMeta

    def __init__(self, inputvis="", cell="", robust=2.0, field="", spw="", stokes="I", M=512, N=512, savemodel=True, verbose=True, **kwargs):
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
                 sidelobethreshold=2.0, minbeamfrac=0.3, deconvolver="hogbom", uvtaper=[], scales=[], uvrange="", pbcor=False, cycleniter=-1, savemodel=True, clean_savemodel=None, **kwargs):
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
        tclean(vis=self.inputvis, imagename=imagename, field=self.field, uvrange=self.uvrange,
               datacolumn=self.datacolumn, specmode=self.specmode, stokes=self.stokes, deconvolver=self.deconvolver, scales=self.scales, nterms=self.nterms,
               imsize=imsize, cell=self.cell, weighting="briggs", robust=self.robust, niter=self.niter, threshold=self.threshold, nsigma=self.nsigma,
               interactive=self.interactive, gridder=self.gridder, pbcor=self.pbcor, uvtaper=self.uvtaper, savemodel=self.clean_savemodel, usemask=self.usemask,
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

    def run(self, imagename=""):
        return


class GPUvmem(Imager):
    def __init(self, executable="gpuvmem", gpublocks=[], initial_values=[], regularization_factors=[], gpu_ids=[], model_in="mod_in.fits", model_out="mod_out.fits", residual_out="residuals.ms", positivity=True, gridding=False, print_images=False, **kwargs):
        super(GPUvmem, self).__unut__(**kwargs)
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(**kwargs)

        make_canvas(model_in)
    def restore(self, restored_image="restored"):
        residual_image = "residual"
        os.system("rm -rf *.log *.last " + residual_image +
                  ".* mod_out convolved_mod_out convolved_mod_out.fits " + restored_image + " " + restored_image + ".fits")

        importfits(imagename="model_out", fitsimage=self.model_out)
        shape = imhead(imagename="model_out", mode="get", hdkey="shape")
        pix_num = shape[0]
        cdelt = imhead(imagename="model_out", mode="get", hdkey="cdelt2")
        cdelta = qa.convert(v=cdelt, outunit="arcsec")
        cdeltd = qa.convert(v=cdelt, outunit="deg")
        pix_size = str(cdelta['value']) + "arcsec"

        tclean(vis=self.residual_out, imagename=residual_image, deconvolver='hogbom', niter=0,
                stokes=self.stokes, weighting='briggs', nterms=1, robust=self.robust, imsize=[self.M, self.N], cell=self.cell datacolumn='RESIDUAL')

        exportfits(imagename=residual_image + ".image",
                   fitsimage=residual_image + ".image.fits", overwrite=True, history=False)

        ia.open(infile=residual_image + ".image")
        rbeam = ia.restoringbeam()
        ia.done()
        ia.close()

        bmaj = imhead(imagename=residual_image + ".image",
                      mode="get", hdkey="beammajor")
        bmin = imhead(imagename=residual_image + ".image",
                      mode="get", hdkey="beamminor")
        bpa = imhead(imagename=residual_image + ".image",
                     mode="get", hdkey="beampa")

        minor = qa.convert(v=bmin, outunit="deg")
        pa = qa.convert(v=bpa, outunit="deg")

        ia.open(infile="model_out")
        ia.convolve2d(outfile="convolved_model_out", axes=[
                      0, 1], type='gauss', major=bmaj, minor=bmin, pa=bpa)
        ia.done()
        ia.close()

        exportfits(imagename="convolved_model_out",
                   fitsimage="convolved_model_out.fits", overwrite=True, history=False)
        ia.open(infile="convolved_model_out.fits")
        ia.setrestoringbeam(beam=rbeam)
        ia.done()
        ia.close()

        imagearr = ["convolved_model_out.fits", residual_image + ".image.fits"]

        immath(imagename=imagearr, expr=" (IM0   + IM1) ", outfile=restored_image)

        exportfits(imagename=restored_image, fitsimage=restored_image +
                   ".fits", overwrite=True, history=False)

        return residual_image + ".image.fits", restored_image + ".fits"

    def make_canvas(self, name="model_input.fits"):
        tclean(vis=self.inputvis, imagename='model_input', specmode='mfs', niter=0,
               deconvolver='hogbom', interactive=False, cell=self.cell, stokes=self.stokes, robust=0.0,
               imsize=[self.M, self.N], weighting='briggs')
        exportfits(imagename='model_input.image', fitsimage=name, overwrite=True);

    def run(self, imagename=""):

        return
