"""
The imager contract of self-calibration (#23). An imager may be any reconstruction method, also a
user's own, as long as it:

1. images the calibrated visibilities, CORRECTED_DATA;
2. writes the model visibilities of its reconstruction to MODEL_DATA of the measurement set it was
   given, because the next gaincal solves against that column.

Imager.image() enforces the second point around every call: it removes any model left from before,
runs the imager, and checks that a model was written for the selected data. The first point can't
be observed cheaply, so snow.testing.ImagerContractSuite proves it on simulated data where only
CORRECTED_DATA is clean.

Plug in an imager by subclassing Imager and implementing _reconstruct, by wrapping a function with
CallableImager, or by name through the "snow.imagers" entry point (snow.imaging.registry).
"""
from __future__ import annotations

import dataclasses
import os
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Callable, Mapping, Optional

from casatasks import delmod

from ..config import Selection
from ..utils.image_utils import read_image
from ..utils.selfcal_utils import is_column_in_ms, model_column_problem
from .statistics import ImageStatistics


class ImagerContractError(RuntimeError):
    """
    An imager broke the self-calibration contract
    """


@dataclass(frozen=True)
class ImagingRequest:
    """
    What to image. The imager receives the measurement set on every call and doesn't store it.

    Parameters
    ----------
    ms :
        Measurement set to image; its MODEL_DATA must hold the reconstruction's model afterwards
    name :
        Prefix of the output images
    selection :
        The data selection shared by the solve, the apply and the imager
    """
    ms: str
    name: str
    selection: Selection = Selection()
    # The contract: imagers image the calibrated data
    data_column: str = field(default="CORRECTED_DATA", init=False)


@dataclass(frozen=True)
class ImagingResult:
    """
    The images of a reconstruction: CASA images (directories) or FITS files

    Parameters
    ----------
    restored :
        Restored image, in Jy/beam
    residual :
        Residual image, in Jy/beam
    model :
        Model image, in Jy/pixel
    mask :
        Image that is non-zero where there is emission, left out of the noise estimate. A mask
        that covers the whole image is ignored
    """
    restored: str
    residual: str
    model: str
    mask: Optional[str] = None

    def statistics(self) -> ImageStatistics:
        """
        The image statistics, computed the same way whatever imager produced the images
        """
        restored, beam = read_image(self.restored)
        residual, _ = read_image(self.residual)
        model, _ = read_image(self.model)
        source_mask = None
        if self.mask is not None and os.path.exists(self.mask):
            source_mask = read_image(self.mask)[0] != 0
            # A mask over the whole image (e.g. clean everywhere) says nothing about the source
            if source_mask.all():
                source_mask = None
        return ImageStatistics.from_arrays(restored, residual, model, source_mask, beam)


class Imager(ABC):
    """
    Base class of the imagers. Subclasses implement _reconstruct; image() is final, so the contract
    is always enforced.
    """

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)
        if "image" in cls.__dict__:
            raise TypeError(
                "{0} overrides Imager.image(), which enforces the self-calibration contract; "
                "implement _reconstruct() instead".format(cls.__name__)
            )

    @abstractmethod
    def _reconstruct(self, request: ImagingRequest) -> ImagingResult:
        """
        Images request.data_column of request.ms with request.selection and writes the model
        visibilities to MODEL_DATA of request.ms
        """

    def with_overrides(self, **overrides) -> Imager:
        """
        A copy of the imager with some options changed, e.g. for one self-calibration step.
        Dataclass imagers get it from dataclasses.replace, which rejects unknown names.
        """
        if not dataclasses.is_dataclass(self):
            raise NotImplementedError(
                "{0} is not a dataclass: implement with_overrides()".format(type(self).__name__)
            )
        return dataclasses.replace(self, **overrides)

    def image(self, request: ImagingRequest) -> ImagingResult:
        """
        Runs the imager and enforces the self-calibration contract

        Raises
        ------
        ImagerContractError
            If the measurement set has no CORRECTED_DATA, or the imager didn't write MODEL_DATA
        """
        if not is_column_in_ms(request.ms, request.data_column):
            raise ImagerContractError(
                "{0} has no {1}: self-calibration images the calibrated data".format(
                    request.ms, request.data_column
                )
            )
        # A model left from an earlier round, real or virtual, can't pass the check below
        delmod(vis=request.ms, otf=True, scr=True)
        result = self._reconstruct(request)
        problem = model_column_problem(request.ms, request.selection.field, request.selection.spw)
        if problem:
            raise ImagerContractError(
                "{0} broke the self-calibration contract: {1}. An imager must write the model "
                "visibilities of its reconstruction to MODEL_DATA of {2}.".format(
                    type(self).__name__, problem, request.ms
                )
            )
        return result


@dataclass(frozen=True)
class CallableImager(Imager):
    """
    Wraps a function as an imager: function(request, **options) images the request and returns an
    ImagingResult. The contract is enforced as for any imager.

    Parameters
    ----------
    function :
        The imaging function
    options :
        Keyword arguments passed to the function; with_overrides() changes them
    """
    function: Callable[..., ImagingResult]
    options: Mapping = field(default_factory=dict, hash=False)

    def __post_init__(self):
        if not callable(self.function):
            raise TypeError("function must be callable, not {0!r}".format(self.function))
        object.__setattr__(self, "options", MappingProxyType(dict(self.options)))

    def with_overrides(self, **overrides) -> CallableImager:
        return dataclasses.replace(self, options={**self.options, **overrides})

    def _reconstruct(self, request: ImagingRequest) -> ImagingResult:
        return self.function(request, **self.options)
