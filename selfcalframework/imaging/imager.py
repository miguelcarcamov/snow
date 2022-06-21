from abc import ABCMeta, abstractmethod

from selfcalframework.utils.image_utils import *


class Imager(metaclass=ABCMeta):

    def __init__(
        self,
        inputvis="",
        output="",
        cell="",
        robust=2.0,
        weighting="briggs",
        field="",
        spw="",
        stokes="I",
        phasecenter="",
        datacolumn="corrected",
        M=512,
        N=512,
        niter=100,
        noise_pixels=None,
        savemodel=True,
        verbose=True
    ):
        self.psnr = 0.0
        self.peak = 0.0
        self.stdv = 0.0
        self.name = ""
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])

        if self.noise_pixels is None:
            self.noise_pixels = -1  # Get to the last pixel for x and y axes
        self.nantennas = calculate_number_antennas(inputvis)
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

    def setVerbose(self, verbose=True):
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

    def calculateStatistics_FITS(
        self, signal_fits_name="", residual_fits_name="", stdv_pixels=None
    ):
        if stdv_pixels is None:
            psnr, peak, stdv = calculatePSNR_FITS(
                signal_fits_name, residual_fits_name, self.noise_pixels
            )
        else:
            psnr, peak, stdv = calculatePSNR_FITS(signal_fits_name, residual_fits_name, stdv_pixels)

        self.psnr = peak / stdv
        self.peak = peak
        self.stdv = stdv

    def calculateStatistics_MSImage(self, signal_ms_name="", residual_ms_name="", stdv_pixels=None):
        if stdv_pixels is None:
            psnr, peak, stdv = calculatePSNR_MS(signal_ms_name, residual_ms_name, self.noise_pixels)

        self.psnr = peak / stdv
        self.peak = peak
        self.stdv = stdv

    @abstractmethod
    def run(self, imagename=""):
        return
