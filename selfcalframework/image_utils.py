from astropy.io import fits
from split import split
from exportfits import exportfits
import numpy as np


def get_data(fits_name=""):
    image_data = fits.getdata(fits_name)
    image_data = np.nan_to_num(image_data)
    return image_data.squeeze()


def exportMStoFITS(msname=""):
    fitsfile_name = msname + ".fits"
    exportfits(msname, fitsfile_name)
    return fitsfile_name


def calculatePSNR_FITS(signal_fits_name="", residual_fits_name="", pixels=80):
    signal_data = get_data(signal_fits_name)
    res_data = get_data(residual_fits_name)
    peak = np.max(np.max(signal_data))
    stdv = np.std(res_data[0:pixels, 0:pixels])
    return peak / stdv, peak, stdv


def calculatePSNR_MS(signal_ms_name="", residual_ms_name="", pixels=80):
    fits_signal = exportMStoFITS(msname=signal_ms_name)
    fits_residual = exportMStoFITS(msname=residual_ms_name)
    psnr, peak, stdv = calculatePSNR_FITS(fits_signal, fits_residual, pixels)
    return psnr, peak, stdv
