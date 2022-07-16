import os
from typing import Tuple, Union
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from casatasks import exportfits, imstat
from reproject import reproject_interp


def nanrms(x, axis=None) -> Union[np.ndarray, float]:
    """
    Function that calculates the root-mean-squared of a numpy array discarding nan values.

    Parameters
    ----------
    x : Input numpy array
    axis : Axis or axes along which the calculation is performed. The default, axis=None, will calculate the RMS using
    all the elements of the input array. If axis is negative it counts from the last to the first axis.

    Returns
    -------
    The RMS of the input numpy array discarding nan values
    """
    return np.sqrt(np.nanmean(x**2, axis=axis))


def rms(x, axis=None) -> Union[np.ndarray, float]:
    """
    Function that calculates the root-mean-squared of a numpy array.

    Parameters
    ----------
    x : Input numpy array
    axis : Axis or axes along which the calculation is performed. The default, axis=None, will calculate the RMS using
    all the elements of the input array. If axis is negative it counts from the last to the first axis.

    Returns
    -------
    The RMS of the input numpy array
    """
    return np.sqrt(np.mean(x**2, axis=axis))


def get_header(fits_name: str = "") -> dict:
    """
    Function that returns the FITS header from a FITS file.

    Parameters
    ----------
    fits_name : Absolute path to the FITS file

    Returns
    -------
    The FITS header as a dictionary
    """
    header = fits.getheader(fits_name)
    return header


def get_hdu(fits_name: str = "") -> fits.PrimaryHDU:
    """
    Function that returns the Primary HDU at position 0 from a FITS file

    Parameters
    ----------
    fits_name : Absolute path to the FITS file

    Returns
    -------
    The Primary HDU of the fits file
    """
    with fits.open(fits_name) as hdul:
        hdu = hdul[0]
    return hdu


def get_hdul(fits_name: str = "") -> fits.HDUList:
    """
    Function that returns the HDU List from a FITS file

    Parameters
    ----------
    fits_name : Absolute path to the FITS file

    Returns
    -------
    The HDU List
    """
    hdul = fits.open(fits_name)
    return hdul


def get_data(fits_name: str = "") -> np.ndarray:
    """
    Function that gets the data as a numpy array from a FITS file

    Parameters
    ----------
    fits_name : Absolute path to the FITS file

    Returns
    -------
    A numpy array with the FITS file data
    """
    image_data = fits.getdata(fits_name)
    return image_data.squeeze()


def get_header_and_data(fits_name: str = "") -> Tuple[dict, np.ndarray]:
    """
    Function that gets the header and the data from a FITS file

    Parameters
    ----------
    fits_name : Absolute path to the FITS file

    Returns
    -------
    A tuple with the header and the data
    """
    with fits.open(fits_name) as hdul:
        header = hdul[0].header
        data = hdul[0].data
    return header, data


def export_ms_to_fits(msname: str = "") -> str:
    """
    Function that export a CASA image file to a FITS image

    Parameters
    ----------
    msname : Absolute path to the CASA image file

    Returns
    -------
    The absolute path to the FITS file
    """
    fitsfile_name = msname + ".fits"
    exportfits(msname, fitsfile_name)
    return fitsfile_name


def calculate_psnr_fits(
    signal_fits_name: str = "",
    residual_fits_name: str = "",
    pixels: int = None
) -> Tuple[float, float, float]:
    """
    Function that calculates the peak signal-to-noise ratio of a reconstruction with images resulting in FITS files.
    The peak is calculated from the restored image, the RMS is calculated in area of the residual image.

    Parameters
    ----------
    signal_fits_name : The absolute path to the restored image
    residual_fits_name : The absolute path to the residual image
    pixels : Number of pixels on where to calculate the RMS

    Returns
    -------
    A tuple with the peak signal-to-noise, the peak and the RMS
    """
    signal_data = get_data(signal_fits_name)
    res_data = get_data(residual_fits_name)

    stdv = nanrms(res_data[0:pixels, 0:pixels])
    peak = np.nanmax(signal_data)
    psnr = peak / stdv

    return psnr, peak, stdv


def calculate_psnr_ms(signal_ms_name: str = "",
                      residual_ms_name: str = "",
                      pixels: int = None) -> Tuple[float, float, float]:
    """
    Function that calculates the peak signal-to-noise ratio of a reconstruction with images resulting in CASA files.
    The peak is calculated from the restored image, the RMS is calculated in area of the residual image.

    Parameters
    ----------
    signal_fits_name : The absolute path to the restored image
    residual_fits_name : The absolute path to the residual image
    pixels : Number of pixels on where to calculate the RMS

    Returns
    -------
    A tuple with the peak signal-to-noise, the peak and the RMS
    """
    box = ""
    if pixels is not None:
        box = "0,0," + str(pixels - 1) + "," + str(pixels - 1)

    stats_signal = imstat(signal_ms_name, box=box)
    stats_residuals = imstat(residual_ms_name, box=box)
    peak = stats_signal["max"][0]
    stdv = stats_residuals["rms"][0]
    psnr = peak / stdv
    return psnr, peak, stdv


def reproject(fits_file_to_resamp: str = "",
              fits_file_model: str = "",
              order: str = "bilinear") -> Union[str, None]:
    """
    Function that reprojects an image if two images does not have the same size or pixel-size

    Parameters
    ----------
    fits_file_to_resamp : Absolute path to the FITS file image that you want to resample
    fits_file_model : The model FITS file image
    order : The order of the interpolation when reprojecting. This can be any of the following strings:
            ‘nearest-neighbor’
            ‘bilinear’
            ‘biquadratic’
            ‘bicubic’
    Returns
    -------
    A string with the absolute path of the reproject image file
    """
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

        same_astrometry_cond = np.fabs((model_dy / mask_dy) - 1.) < 1E-3

        if same_astrometry_cond or mask_M != model_M or mask_N != model_N:
            print("The mask header is not the same as the model image, resampling...")
            reprojected_array = reproject_interp(
                (data_mask, mask_WCS),
                model_WCS,
                return_footprint=False,
                order=order,
                shape_out=(model_M, model_N)
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
