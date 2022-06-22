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
        phase_center="",
        data_column="corrected",
        M=512,
        N=512,
        niter=100,
        noise_pixels=None,
        save_model=True,
        verbose=True
    ):
        self.inputvis = inputvis
        self.output = output
        self.cell = cell
        self.robust = robust
        self.weighting = weighting
        self.field = field
        self.spw = spw
        self.stokes = stokes
        self.phase_center = phase_center
        self.data_column = data_column
        self.M = M
        self.N = N
        self.niter = niter
        self.noise_pixels = noise_pixels
        self.save_model = save_model
        self.verbose = verbose
        self.psnr = 0.0
        self.peak = 0.0
        self.stdv = 0.0
        self.name = ""

        if self.noise_pixels is None:
            self.noise_pixels = -1  # Get to the last pixel for x and y axes
        self.nantennas = calculate_number_antennas(inputvis)
        # self.__dict__.update(kwargs)
    def _calculateStatistics_FITS(
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

    def _calculateStatistics_MSImage(
        self, signal_ms_name="", residual_ms_name="", stdv_pixels=None
    ):
        if stdv_pixels is None:
            psnr, peak, stdv = calculatePSNR_MS(signal_ms_name, residual_ms_name, self.noise_pixels)

        self.psnr = peak / stdv
        self.peak = peak
        self.stdv = stdv

    @abstractmethod
    def run(self, imagename=""):
        return
