from astropy.io import fits
from casatasks import exportfits
import numpy as np
from casatools import table
import os

tb = table()

def get_data(fits_name=""):
    image_data = fits.getdata(fits_name)
    image_data = np.nan_to_num(image_data)
    return image_data.squeeze()


def exportMStoFITS(msname=""):
    fitsfile_name = msname + ".fits"
    exportfits(msname, fitsfile_name)
    return fitsfile_name


def calculatePSNR_FITS(signal_fits_name="", residual_fits_name="", pixels=100):
    signal_data = get_data(signal_fits_name)
    res_data = get_data(residual_fits_name)
    peak = np.max(signal_data)
    stdv = np.std(res_data[0:pixels, 0:pixels])
    return peak / stdv, peak, stdv


def calculatePSNR_MS(signal_ms_name="", residual_ms_name="", pixels=100):
    fits_signal = exportMStoFITS(msname=signal_ms_name)
    fits_residual = exportMStoFITS(msname=residual_ms_name)
    psnr, peak, stdv = calculatePSNR_FITS(fits_signal, fits_residual, pixels)
    return psnr, peak, stdv

def calculate_number_antennas(ms_name=""):
    if ms_name != "":
        if os.path.exists(ms_name):
            tb.open(tablename=ms_name+"/ANTENNA")
            query_table = tb.taql("select NAME from "+ms_name+"/ANTENNA"+" where !FLAG_ROW")
            nrows = len(query_table.getcol("NAME"))
            tb.close()
            return nrows
        else:
            raise FileNotFoundError("The Measurement Set File does not exist")
    else:
        raise ValueError("Measurement Set File cannot be empty")
