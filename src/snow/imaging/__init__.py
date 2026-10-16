"""
Imagers for self-calibration (#23). Any imager that follows the contract of
snow.imaging.contract can be used: a subclass of Imager, a function wrapped in CallableImager, or
an imager registered by name under the "snow.imagers" entry point.
"""
from .backends import Tclean
from .contract import (CallableImager, Imager, ImagerContractError, ImagingRequest, ImagingResult)
from .registry import available_imagers, imager_from_name
from .statistics import ImageStatistics
