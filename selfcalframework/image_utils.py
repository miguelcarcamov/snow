from astropy.io import fits
from split import split
from exportfits import exportfits
import numpy as np

def get_data(fits_name=""):
    image_data = fits.getdata(name)
    image_data = np.nan_to_num(image_data)
    return image_data.squeeze()

def exportMStoFITS(msname=""):
    fitsfile_name = msname + ".fits"
    exportfits(msname, fits_name)

def calculatePSNR_FITS(signal_fits_name="", residual_fits_name="", pixels=80):
    signal_data = get_data(signal_fits_name)
    res_data = get_data(residual_fits_name)
    peak = np.max(max(data))
    stdv = np.std(res_data[0:pixels,0:pixels])
    return peak/stdv, peak, stdv
