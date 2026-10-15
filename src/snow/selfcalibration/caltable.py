from dataclasses import dataclass
from typing import Sequence


@dataclass(frozen=True)
class CalTable:
    """
    Calibration table accepted by a self-calibration iteration

    Parameters
    ----------
    path :
        Path to the CASA calibration table
    calmode :
        Calibration mode the table was solved with: "p", "a" or "ap"
    spwmap :
        Spectral window map used to apply the table. Empty when every spectral window has its own
        solutions
    """
    path: str
    calmode: str
    spwmap: tuple = ()

    @property
    def solves_phase(self) -> bool:
        return "p" in self.calmode

    @property
    def solves_amplitude(self) -> bool:
        return "a" in self.calmode


def derive_spwmap(combine: str, spw_ids: Sequence[int], nspw: int) -> list:
    """
    Function that derives the spectral window map needed to apply a gain solution. When the
    solution combines spectral windows, CASA stores it under the lowest selected spectral window,
    so every spectral window maps to it. Otherwise each spectral window has its own solution and
    no map is needed.

    Parameters
    ----------
    combine :
        Data axes combined when solving, e.g. "spw" or "scan,spw"
    spw_ids :
        Spectral window ids selected when solving
    nspw :
        Total number of spectral windows in the measurement set

    Returns
    -------
    The spectral window map, or an empty list if the solution does not combine spectral windows
    """
    if "spw" not in [axis.strip() for axis in combine.split(",")]:
        return []
    if not spw_ids:
        raise ValueError("No spectral windows selected to solve for")
    return [min(spw_ids)] * nspw


def build_gaintable_arguments(caltables: Sequence[CalTable], interp: str) -> dict:
    """
    Function that builds the gaintable, spwmap and interp arguments of gaincal (tables pre-applied
    while solving) and applycal (tables applied to the data), one entry per table

    Parameters
    ----------
    caltables :
        Calibration tables in the order they are applied
    interp :
        Temporal interpolation used for every table

    Returns
    -------
    Dictionary with the gaintable, spwmap and interp keyword arguments
    """
    return {
        "gaintable": [caltable.path for caltable in caltables],
        "spwmap": [list(caltable.spwmap) for caltable in caltables],
        "interp": [interp] * len(caltables),
    }
