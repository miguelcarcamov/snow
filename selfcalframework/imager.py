import numpy as np
import os
import time
from casatasks import tclean
from casatasks import fixvis
from casatasks import importfits
from casatasks import exportfits
from casatasks import imhead
from casatasks import immath
from .image_utils import *
from casatools import image
from casatools import quanta
from abc import ABCMeta, abstractmethod
import shlex
import subprocess

class Imager(metaclass=ABCMeta):

    def __init__(self, inputvis="", output="", cell="", robust=2.0, weighting="briggs", field="", spw="", stokes="I", phasecenter="", datacolumn="corrected", M=512, N=512, niter=100, savemodel=True, verbose=True):
        self.psnr = 0.0
        self.peak = 0.0
        self.stdv = 0.0
        self.name = ""
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        # self.__dict__.update(kwargs)

    def getVis(self):
        return self.inputvis

    def setVis(self, inputvis=""):
        self.inputvis = inputvis

    def getCell(self):
        return self.cell

    def setCell(self, cell=""):
        self.cell = cell

    def getRobust(self):
        return self.robust

    def setRobust(self, robust=0.0):
        self.robust = robust

    def getField(self):
        return self.field

    def setField(self, field=""):
        self.field = field

    def getSpw(self):
        return self.spw

    def setSpw(self, spw=""):
        self.spw = spw

    def getStokes(self):
        return self.stokes

    def setStokes(self, stokes=""):
        self.stokes = stokes

    def getMN(self):
        return self.M, self.N

    def setMN(self, M, N):
        self.M = M
        self.N = N

    def getSaveModel(self):
        return self.savemodel

    def setSaveModel(self, savemodel=False):
        self.savemodel = savemodel

    def getVerbose(self):
        return self.verbose

    def setVerbose(self):
        self.verbose = verbose

    def getOutputPath(self):
        return self.output

    def getPSNR(self):
        return self.psnr

    def getPeak(self):
        return self.peak

    def getSTDV(self):
        return self.stdv

    def getPhaseCenter(self):
        return self.phasecenter

    def setPhaseCenter(self, phasecenter=""):
        self.phasecenter = phasecenter

    def calculateStatistics_FITS(self, signal_fits_name="", residual_fits_name="", stdv_pixels=80):
        self.psnr, self.peak, self.stdv = calculatePSNR_FITS(
            signal_fits_name, residual_fits_name, stdv_pixels)

    def calculateStatistics_MSImage(self, signal_ms_name="", residual_ms_name="", stdv_pixels=80):
        self.psnr, self.peak, self.stdv = calculatePSNR_MS(
            signal_ms_name, residual_ms_name, stdv_pixels)

    @abstractmethod
    def run(self, imagename=""):
        return


class Clean(Imager):
    def __init__(self, nterms=1, threshold=0.0, nsigma=0.0, interactive=False, mask="", usemask="auto-multithresh", negativethreshold=0.0, lownoisethreshold=1.5, noisethreshold=4.25,
                 sidelobethreshold=2.0, minbeamfrac=0.3, specmode="", gridder="standard", wprojplanes=-1, deconvolver="hogbom", uvtaper=[], scales=[], uvrange="", pbcor=False, cycleniter=0, clean_savemodel=None, **kwargs):
        super(Clean, self).__init__(**kwargs)
        self.name = "TClean"
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        # self.__dict__.update(kwargs)

        if(self.savemodel):
            self.clean_savemodel = "modelcolumn"

    def getNSigma(self):
        return self.nsigma
    def setNSigma(self, nsigma):
        self.nsigma = nsigma

    def getThreshold(self):
        return self.threshold
    def setThreshold(self):
        return self.threshold

    def run(self, imagename=""):
        imsize = [self.M, self.N]
        tclean(vis=self.inputvis, imagename=imagename, field=self.field, phasecenter=self.phasecenter, uvrange=self.uvrange,
               datacolumn=self.datacolumn, specmode=self.specmode, stokes=self.stokes, deconvolver=self.deconvolver, scales=self.scales, nterms=self.nterms,
               imsize=imsize, cell=self.cell, weighting=self.weighting, robust=self.robust, niter=self.niter, threshold=self.threshold, nsigma=self.nsigma,
               interactive=self.interactive, gridder=self.gridder, mask=self.mask, pbcor=self.pbcor, uvtaper=self.uvtaper, savemodel=self.clean_savemodel, usemask=self.usemask,
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


class GPUvmem(Imager):
    def __init__(self, executable="gpuvmem", gpublocks=[16, 16, 256], initialvalues=[], regfactors=[], gpuids=[0], residualoutput="residuals.ms",
                 modelin="mod_in.fits", modelout="mod_out.fits", griddingthreads=4, positivity=True, ftol=1e-12, noise_cut = 10.0, gridding=False, printimages=False, **kwargs):
        super(GPUvmem, self).__init__(**kwargs)
        self.name = "GPUvmem"
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        # self.__dict__.update(kwargs)
        if(self.phasecenter != ""):
            fixvis(vis=self.visfile, outputvis=self.visfile, field=self.field, phasecenter=self.phasecenter)

    def getRegfactors(self):
        return self.regfactors

    def setRegFactors(self, regfactors):
        self.factors = regfactors

    def _restore(self, model_fits="", residual_ms="", restored_image="restored"):
        qa = quanta()
        ia = image()
        residual_image = "residual"
        os.system("rm -rf *.log *.last " + residual_image +
                  ".* mod_out convolved_mod_out convolved_mod_out.fits " + restored_image + " " + restored_image + ".fits")

        importfits(imagename="model_out", fitsimage=model_fits, overwrite=True)
        shape = imhead(imagename="model_out", mode="get", hdkey="shape")
        pix_num = shape[0]
        cdelt = imhead(imagename="model_out", mode="get", hdkey="cdelt2")
        cdelta = qa.convert(v=cdelt, outunit="arcsec")
        cdeltd = qa.convert(v=cdelt, outunit="deg")
        pix_size = str(cdelta['value']) + "arcsec"

        tclean(vis=residual_ms, imagename=residual_image, specmode='mfs', deconvolver='hogbom', niter=0,
               stokes=self.stokes, nterms=1, weighting=self.weighting, robust=self.robust, imsize=[self.M, self.N], cell=self.cell, datacolumn='data')

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
        im2 = ia.convolve2d(outfile="convolved_model_out", axes=[
                      0, 1], type='gauss', major=bmaj, minor=bmin, pa=bpa, overwrite=True)
        im2.done()
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

    def _make_canvas(self, name="model_input"):
        fitsimage = name + '.fits'
        tclean(vis=self.inputvis, imagename=name, specmode='mfs', niter=0,
               deconvolver='hogbom', interactive=False, cell=self.cell, stokes=self.stokes, robust=self.robust,
               imsize=[self.M, self.N], weighting=self.weighting)
        exportfits(imagename=name + '.image',
                   fitsimage=fitsimage, overwrite=True)
        return fitsimage

    def run(self, imagename=""):
        model_input = self._make_canvas(imagename + "_input")
        model_output = imagename + ".fits"
        residual_output = imagename + "_" + self.residualoutput
        restored_image = imagename + ".restored"

        args = self.executable + " -X " + str(self.gpublocks[0]) + " -Y " + str(self.gpublocks[1]) + " -V " + str(self.gpublocks[2]) \
            + " -i " + self.inputvis + " -o " + residual_output + " -z " + ",".join(map(str, self.initialvalues)) \
            + " -Z " + ",".join(map(str, self.regfactors)) + " -G " + ",".join(map(str, self.gpuids)) \
            + " -m " + model_input + " -O " + model_output + " -N " + str(self.noise_cut) \
            + " -R " + str(self.robust) + " -t " + str(self.niter)

        if(self.gridding):
            args += " -g " + str(self.griddingthreads)

        if(self.printimages):
            args += " --print-images"

        if(not self.positivity):
            args += " --nopositivity"

        if(self.verbose):
            args += " --verbose"

        if(self.savemodel):
            args += " --savemodel-input"

        print(args)
        args = shlex.split(args)
        print(args)

        # Run gpuvmem and wait until it finishes
        p = subprocess.Popen(args, env=os.environ)
        p.wait()

        # Restore the image
        residual_fits, restored_fits = self._restore(model_fits=model_output,
                                                     residual_ms=residual_output, restored_image=restored_image)

        # Calculate SNR and standard deviation
        self.calculateStatistics_FITS(
            signal_fits_name=restored_fits, residual_fits_name=residual_fits)


class WSClean(Imager):
    def __init__(self, **kwargs):
        super(WSClean, self).__init__(**kwargs)
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(kwargs)

    def run(self, imagename=""):
        return
