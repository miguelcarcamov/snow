"""
CASA tclean as a self-calibration imager (#23). Options pass through with tclean's own names and are
checked against its signature; the ones that self-calibration controls come from the request, the
image geometry and the weighting, so each is set in exactly one place.
"""
from __future__ import annotations

import dataclasses
import inspect
import os
from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Mapping

from casatasks import tclean

from ...config import ImageGeometry, Weighting
from ..contract import Imager, ImagingRequest, ImagingResult

# tclean parameters set by snow, and where each comes from
CONTROLLED = {
    "vis": "the request",
    "imagename": "the request",
    "field": "the Selection",
    "spw": "the Selection",
    "uvrange": "the Selection",
    "datacolumn": "the contract (corrected)",
    "savemodel": "the contract (modelcolumn)",
    "imsize": "ImageGeometry",
    "cell": "ImageGeometry",
    "phasecenter": "ImageGeometry",
    "stokes": "ImageGeometry",
    "weighting": "Weighting",
    "robust": "Weighting",
    "uvtaper": "Weighting",
}


def _tclean_parameters() -> frozenset:
    return frozenset(inspect.signature(tclean).parameters)


@dataclass(frozen=True, init=False)
class Tclean(Imager):
    """
    Images with CASA tclean, e.g.

        Tclean(ImageGeometry(imsize=1024, cell="0.05arcsec"), Weighting("briggs", robust=0.5),
               deconvolver="hogbom", niter=1000, usemask="auto-multithresh", nsigma=3.0)

    Parameters
    ----------
    geometry :
        Image grid; imsize and cell are required
    weighting :
        Visibility weighting
    tclean_options :
        Any other tclean parameter, with tclean's name
    """
    geometry: ImageGeometry
    weighting: Weighting
    options: Mapping = field(hash=False)

    def __init__(
        self, geometry: ImageGeometry, weighting: Weighting = Weighting(), **tclean_options
    ):
        if not isinstance(geometry, ImageGeometry) or not isinstance(weighting, Weighting):
            raise TypeError("Tclean needs an ImageGeometry and a Weighting")
        if geometry.imsize is None or not geometry.cell:
            raise ValueError(
                "Tclean needs the image size and the cell: ImageGeometry(imsize=..., cell=...)"
            )
        _check_options(tclean_options)
        object.__setattr__(self, "geometry", geometry)
        object.__setattr__(self, "weighting", weighting)
        object.__setattr__(self, "options", MappingProxyType(dict(tclean_options)))

    def with_overrides(self, **overrides) -> Tclean:
        """
        A copy with some options changed: geometry and weighting replace the config objects, any
        other name is a tclean parameter
        """
        geometry = overrides.pop("geometry", self.geometry)
        weighting = overrides.pop("weighting", self.weighting)
        return Tclean(geometry, weighting, **{**self.options, **overrides})

    def _image_suffix(self) -> str:
        # Multi-term deconvolution writes one image per Taylor term; the zeroth is the intensity
        return ".tt0" if self.options.get("deconvolver") == "mtmfs" else ""

    def _reconstruct(self, request: ImagingRequest) -> ImagingResult:
        tclean(
            vis=request.ms,
            imagename=request.name,
            **dataclasses.asdict(request.selection),
            datacolumn="corrected",
            savemodel="modelcolumn",
            imsize=list(self.geometry.imsize),
            cell=self.geometry.cell,
            phasecenter=self.geometry.phasecenter,
            stokes=self.geometry.stokes,
            weighting=self.weighting.weighting,
            robust=self.weighting.robust,
            uvtaper=list(self.weighting.uvtaper),
            **self.options
        )
        suffix = self._image_suffix()
        result = ImagingResult(
            restored=request.name + ".image" + suffix,
            residual=request.name + ".residual" + suffix,
            model=request.name + ".model" + suffix,
            mask=request.name + ".mask",
        )
        for path in (result.restored, result.residual, result.model):
            if not os.path.exists(path):
                raise RuntimeError("tclean did not write {0}".format(path))
        return result


def _check_options(options: Mapping) -> None:
    controlled = sorted(set(options) & set(CONTROLLED))
    if controlled:
        raise TypeError(
            "Tclean sets {0} itself; they come from {1}".format(
                controlled, ", ".join(sorted({CONTROLLED[name]
                                              for name in controlled}))
            )
        )
    unknown = sorted(set(options) - _tclean_parameters())
    if unknown:
        raise TypeError("{0} are not tclean parameters".format(unknown))
