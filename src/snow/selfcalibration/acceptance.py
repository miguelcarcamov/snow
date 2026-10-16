"""
Acceptance policies (#28): whether a self-calibration round is kept or rolled back. A higher PSNR
alone is not enough: amplitude gains can absorb noise and unmodelled flux, which lowers the
residual rms while the image gets worse. The default policy therefore also checks that no model
flux is lost, that the beam is stable and that few solutions were flagged.

Policies decide and never act: the loop (#22c) restores the data when a round is rejected.
"""
from __future__ import annotations

import math
from abc import ABC, abstractmethod
from dataclasses import dataclass
from typing import Tuple

from ..imaging.statistics import ImageStatistics


def _require_in_range(
    policy: object, name: str, low: float, high: float, include_low: bool = True
) -> None:
    value = getattr(policy, name)
    above_low = value >= low if include_low else value > low
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not above_low \
            or not value <= high:
        raise ValueError(
            "{0}.{1} must be in {2}{3}, {4}], not {5!r}".format(
                type(policy).__name__, name, "[" if include_low else "(", low, high, value
            )
        )


@dataclass(frozen=True)
class RoundMetrics:
    """
    What a round is judged on

    Parameters
    ----------
    statistics :
        Statistics of the round's image
    flagged_fraction :
        Fraction of the round's gain solutions that were flagged (0 for the initial image)
    """
    statistics: ImageStatistics
    flagged_fraction: float = 0.0

    def __post_init__(self):
        if not 0.0 <= self.flagged_fraction <= 1.0:
            raise ValueError(
                "flagged_fraction must be between 0 and 1, not {0!r}".format(self.flagged_fraction)
            )


@dataclass(frozen=True)
class Decision:
    """
    Whether to keep a round, and why
    """
    accepted: bool
    reasons: Tuple[str, ...] = ()

    @classmethod
    def accept(cls, reason: str) -> Decision:
        return cls(True, (reason, ))

    @classmethod
    def reject(cls, reason: str) -> Decision:
        return cls(False, (reason, ))


class AcceptancePolicy(ABC):
    """
    Decides whether the current round is kept, comparing it with the last accepted one
    """

    @abstractmethod
    def decide(self, previous: RoundMetrics, current: RoundMetrics) -> Decision:
        ...


@dataclass(frozen=True)
class PSNRImproves(AcceptancePolicy):
    """
    The peak signal-to-noise ratio improves by more than a margin

    Parameters
    ----------
    min_ratio :
        The current PSNR must be above min_ratio times the previous one
    """
    min_ratio: float = 1.0

    def __post_init__(self):
        _require_in_range(self, "min_ratio", 0.0, math.inf, include_low=False)

    def decide(self, previous: RoundMetrics, current: RoundMetrics) -> Decision:
        before, after = previous.statistics.psnr, current.statistics.psnr
        message = "PSNR {0:.3f} -> {1:.3f} (needs > {2:.3f})".format(
            before, after, before * self.min_ratio
        )
        if after > before * self.min_ratio:
            return Decision.accept(message)
        return Decision.reject(message)


@dataclass(frozen=True)
class FluxNotLost(AcceptancePolicy):
    """
    The model keeps its flux. Phase self-calibration usually raises it (less decorrelation); a
    drop is the sign of amplitude gains absorbing flux

    Parameters
    ----------
    max_loss :
        Largest allowed fractional drop of the model flux
    """
    max_loss: float = 0.1

    def __post_init__(self):
        _require_in_range(self, "max_loss", 0.0, 1.0)

    def decide(self, previous: RoundMetrics, current: RoundMetrics) -> Decision:
        before, after = previous.statistics.model_flux, current.statistics.model_flux
        message = "model flux {0:.4g} -> {1:.4g} Jy (max loss {2:.0%})".format(
            before, after, self.max_loss
        )
        if after >= (1.0 - self.max_loss) * before:
            return Decision.accept(message)
        return Decision.reject(message)


@dataclass(frozen=True)
class BeamStable(AcceptancePolicy):
    """
    The restoring beam axes change by at most a fraction, e.g. because flagging removed antennas.
    Rounds without beam information pass

    Parameters
    ----------
    tolerance :
        Largest allowed fractional change of the major and minor axes
    """
    tolerance: float = 0.05

    def __post_init__(self):
        _require_in_range(self, "tolerance", 0.0, math.inf)

    def decide(self, previous: RoundMetrics, current: RoundMetrics) -> Decision:
        before, after = previous.statistics.beam, current.statistics.beam
        if before is None or after is None:
            return Decision.accept("no beam to compare")
        change = max(abs(a - b) / b for a, b in zip(after[:2], before[:2]))
        message = "beam axes change {0:.1%} (max {1:.0%})".format(change, self.tolerance)
        if change <= self.tolerance:
            return Decision.accept(message)
        return Decision.reject(message)


@dataclass(frozen=True)
class FlaggedFractionBelow(AcceptancePolicy):
    """
    Few of the round's gain solutions were flagged (e.g. below minsnr)

    Parameters
    ----------
    max_fraction :
        Largest allowed fraction of flagged solutions
    """
    max_fraction: float = 0.25

    def __post_init__(self):
        _require_in_range(self, "max_fraction", 0.0, 1.0)

    def decide(self, previous: RoundMetrics, current: RoundMetrics) -> Decision:
        message = "{0:.1%} of the solutions flagged (max {1:.0%})".format(
            current.flagged_fraction, self.max_fraction
        )
        if current.flagged_fraction <= self.max_fraction:
            return Decision.accept(message)
        return Decision.reject(message)


@dataclass(frozen=True, init=False)
class AllOf(AcceptancePolicy):
    """
    Accepts when every policy accepts. All policies are evaluated, so the decision lists every
    reason
    """
    policies: Tuple[AcceptancePolicy, ...]

    def __init__(self, *policies: AcceptancePolicy):
        if not policies:
            raise ValueError("AllOf needs at least one policy")
        object.__setattr__(self, "policies", policies)

    def decide(self, previous: RoundMetrics, current: RoundMetrics) -> Decision:
        decisions = [policy.decide(previous, current) for policy in self.policies]
        return Decision(
            all(decision.accepted for decision in decisions),
            tuple(reason for decision in decisions for reason in decision.reasons)
        )


@dataclass(frozen=True, init=False)
class AnyOf(AcceptancePolicy):
    """
    Accepts when at least one policy accepts. All policies are evaluated, so the decision lists
    every reason
    """
    policies: Tuple[AcceptancePolicy, ...]

    def __init__(self, *policies: AcceptancePolicy):
        if not policies:
            raise ValueError("AnyOf needs at least one policy")
        object.__setattr__(self, "policies", policies)

    def decide(self, previous: RoundMetrics, current: RoundMetrics) -> Decision:
        decisions = [policy.decide(previous, current) for policy in self.policies]
        return Decision(
            any(decision.accepted for decision in decisions),
            tuple(reason for decision in decisions for reason in decision.reasons)
        )


# PSNR never decides alone (selfcal-domain.mdc)
DEFAULT_ACCEPTANCE = AllOf(PSNRImproves(), FluxNotLost(), BeamStable(), FlaggedFractionBelow())


def solution_snr(
    flux: float,
    sigma_visibility: float,
    nantennas: int,
    nintegrations: int,
    nchannels: int = 1,
    npolarizations: int = 2
) -> float:
    """
    Pre-flight estimate of the signal-to-noise ratio of one antenna's gain solution: each antenna
    is solved from its nantennas - 1 baselines, averaged over the integrations, channels and
    polarizations combined in the solution interval

        S/N = (S / sigma_vis) * sqrt(nantennas - 1) * sqrt(nintegrations * nchannels * npolarizations)

    Compare it with minsnr to check that a solint can be supported before solving.

    Parameters
    ----------
    flux :
        Flux of the model, in Jy
    sigma_visibility :
        Noise of one visibility (one integration, channel and polarization), in Jy
    nantennas :
        Number of antennas
    nintegrations :
        Integrations in one solution interval
    nchannels :
        Channels combined in one solution
    npolarizations :
        Polarizations combined in one solution (2 for gaintype "T", 1 for "G")
    """
    if sigma_visibility <= 0:
        raise ValueError("sigma_visibility must be positive")
    if nantennas < 2 or min(nintegrations, nchannels, npolarizations) < 1:
        raise ValueError("Need at least 2 antennas and 1 integration, channel and polarization")
    return flux / sigma_visibility * math.sqrt(nantennas - 1) * math.sqrt(
        nintegrations * nchannels * npolarizations
    )
