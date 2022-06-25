import os
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from casatasks import exportfits, imstat
from reproject import reproject_interp


def nanrms(x, axis=None):
    return np.sqrt(np.nanmean(x**2, axis=axis))


def rms(x, axis=None):
    return np.sqrt(np.mean(x**2, axis=axis))


def get_header(fits_name: str = ""):
    header = fits.getheader(fits_name)
    return header


def get_hdu(fits_name: str = ""):
    with fits.open(fits_name) as hdul:
        hdu = hdul[0]
    return hdu


def get_hdul(fits_name: str = ""):
    hdul = fits.open(fits_name)
    return hdul


def get_data(fits_name: str = ""):
    image_data = fits.getdata(fits_name)
    return image_data.squeeze()


def get_header_and_data(fits_name: str = ""):
    with fits.open(fits_name) as hdul:
        header = hdul[0].header
        data = hdul[0].data
    return header, data


def export_ms_to_fits(msname: str = ""):
    fitsfile_name = msname + ".fits"
    exportfits(msname, fitsfile_name)
    return fitsfile_name


def calculate_psnr_fits(
    signal_fits_name: str = "", residual_fits_name: str = "", pixels: int = None
):
    signal_data = get_data(signal_fits_name)
    res_data = get_data(residual_fits_name)

    stdv = nanrms(res_data[0:pixels, 0:pixels])
    peak = np.nanmax(signal_data)
    psnr = peak / stdv

    return psnr, peak, stdv


def calculate_psnr_ms(signal_ms_name: str = "", residual_ms_name: str = "", pixels: int = None):
    box = ""
    if pixels is not None:
        box = "0,0," + str(pixels - 1) + "," + str(pixels - 1)

    stats_signal = imstat(signal_ms_name, box=box)
    stats_residuals = imstat(residual_ms_name, box=box)
    peak = stats_signal["max"][0]
    stdv = stats_residuals["rms"][0]
    psnr = peak / stdv
    return psnr, peak, stdv


def reproject(fits_file_to_resamp="", fits_file_model="", order="bilinear"):
    if os.path.exists(fits_file_to_resamp) and os.path.exists(fits_file_model):
        header_mask = get_header(fits_file_to_resamp)
        data_mask = get_data(fits_file_to_resamp)
        header_model = get_header(fits_file_model)
        model_WCS = WCS(header=header_model, naxis=2)
        mask_WCS = WCS(header=header_mask, naxis=2)

        model_M = header_model['NAXIS1']
        model_N = header_model['NAXIS2']
        model_dy = header_model['CDELT2']

        mask_M = header_mask['NAXIS1']
        mask_N = header_mask['NAXIS2']
        mask_dy = header_mask['CDELT2']

        if np.fabs(
            (model_dy - mask_dy) / model_dy
        ) < 1E-3 or mask_M != model_M or mask_N != model_N:
            print("The mask header is not the same as the model image, resampling...")
            reprojected_array, footprint = reproject_interp(
                (data_mask, mask_WCS), model_WCS, order=order, shape_out=(model_M, model_N)
            )
            path_object = Path(fits_file_to_resamp)
            resampled_mask_name = "{0}_{2}{1}".format(
                Path.joinpath(path_object.parent, path_object.stem), path_object.suffix, "resampled"
            )
            fits.writeto(resampled_mask_name, reprojected_array, header_model, overwrite=True)
            return resampled_mask_name
        else:
            return None
    else:
        print("Fits file to reproject: {0}".format(fits_file_to_resamp))
        print("Fits file model: {0}".format(fits_file_model))
        raise FileNotFoundError("Either user mask of model input does not exist")
