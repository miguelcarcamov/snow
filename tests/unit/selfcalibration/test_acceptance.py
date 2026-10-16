"""Acceptance policies (#28): composable checks, and PSNR never decides alone."""
import math

import pytest

from snow.imaging.statistics import ImageStatistics
from snow.selfcalibration.acceptance import (
    DEFAULT_ACCEPTANCE, AllOf, AnyOf, BeamStable, FlaggedFractionBelow, FluxNotLost, PSNRImproves,
    RoundMetrics, solution_snr
)


def metrics(peak=1.0, rms=1e-3, flux=1.0, beam=(1.0, 0.8, 30.0), flagged=0.0) -> RoundMetrics:
    return RoundMetrics(ImageStatistics(peak, rms, flux, beam), flagged)


BEFORE = metrics()


class TestPSNRImproves:

    def test_higher_psnr_is_accepted(self):
        assert PSNRImproves().decide(BEFORE, metrics(rms=0.5e-3)).accepted

    def test_equal_psnr_is_rejected(self):
        assert not PSNRImproves().decide(BEFORE, metrics()).accepted

    def test_margin(self):
        policy = PSNRImproves(min_ratio=1.1)

        assert not policy.decide(BEFORE, metrics(peak=1.05)).accepted
        assert policy.decide(BEFORE, metrics(peak=1.2)).accepted


class TestFluxNotLost:

    def test_more_flux_after_phase_calibration_is_fine(self):
        assert FluxNotLost().decide(BEFORE, metrics(flux=1.3)).accepted

    def test_small_loss_within_tolerance_is_fine(self):
        assert FluxNotLost(max_loss=0.1).decide(BEFORE, metrics(flux=0.95)).accepted

    def test_flux_absorbed_by_the_gains_is_rejected(self):
        decision = FluxNotLost(max_loss=0.1).decide(BEFORE, metrics(flux=0.8))

        assert not decision.accepted and "model flux" in decision.reasons[0]


class TestBeamStable:

    def test_small_change_is_fine(self):
        assert BeamStable(0.05).decide(BEFORE, metrics(beam=(1.03, 0.8, 35.0))).accepted

    def test_large_change_is_rejected(self):
        assert not BeamStable(0.05).decide(BEFORE, metrics(beam=(1.0, 0.9, 30.0))).accepted

    def test_rounds_without_a_beam_pass(self):
        assert BeamStable().decide(BEFORE, metrics(beam=None)).accepted


class TestFlaggedFractionBelow:

    def test_few_flagged_solutions_are_fine(self):
        assert FlaggedFractionBelow(0.25).decide(BEFORE, metrics(flagged=0.1)).accepted

    def test_many_flagged_solutions_are_rejected(self):
        assert not FlaggedFractionBelow(0.25).decide(BEFORE, metrics(flagged=0.4)).accepted


class TestCombinators:

    def test_all_of_needs_every_policy_and_lists_every_reason(self):
        decision = AllOf(PSNRImproves(),
                         FluxNotLost()).decide(BEFORE, metrics(rms=0.5e-3, flux=0.5))

        assert not decision.accepted
        assert len(decision.reasons) == 2

    def test_any_of_needs_one_policy(self):
        policy = AnyOf(PSNRImproves(), FluxNotLost())

        assert policy.decide(BEFORE, metrics(flux=1.2)).accepted
        assert not policy.decide(BEFORE, metrics(rms=2e-3, flux=0.5)).accepted

    def test_combinators_nest_and_are_hashable(self):
        policy = AllOf(PSNRImproves(), AnyOf(FluxNotLost(), BeamStable()))

        assert hash(policy) == hash(AllOf(PSNRImproves(), AnyOf(FluxNotLost(), BeamStable())))

    @pytest.mark.parametrize("combinator", [AllOf, AnyOf])
    def test_combinators_need_a_policy(self, combinator):
        with pytest.raises(ValueError):
            combinator()


class TestDefaultPolicy:

    def test_a_better_image_is_accepted(self):
        assert DEFAULT_ACCEPTANCE.decide(
            BEFORE, metrics(rms=0.5e-3, flux=1.1, flagged=0.05)
        ).accepted

    def test_higher_psnr_with_flux_absorbed_is_rejected(self):
        """The failure PSNR-only acceptance misses: amplitude gains absorbing flux and noise."""
        current = metrics(peak=0.9, rms=0.4e-3, flux=0.7)
        assert current.statistics.psnr > BEFORE.statistics.psnr

        decision = DEFAULT_ACCEPTANCE.decide(BEFORE, current)

        assert not decision.accepted

    def test_higher_psnr_with_most_solutions_flagged_is_rejected(self):
        assert not DEFAULT_ACCEPTANCE.decide(BEFORE, metrics(rms=0.5e-3, flagged=0.6)).accepted


@pytest.mark.parametrize(
    "build", [
        lambda: PSNRImproves(min_ratio=0),
        lambda: FluxNotLost(max_loss=1.5),
        lambda: BeamStable(tolerance=-0.1),
        lambda: FlaggedFractionBelow(max_fraction=2),
        lambda: RoundMetrics(BEFORE.statistics, flagged_fraction=1.5),
    ]
)
def test_invalid_parameters_raise(build):
    with pytest.raises(ValueError):
        build()


class TestSolutionSNR:

    def test_formula(self):
        snr = solution_snr(
            flux=1.0,
            sigma_visibility=0.5,
            nantennas=10,
            nintegrations=6,
            nchannels=4,
            npolarizations=2
        )

        assert snr == pytest.approx(2.0 * 3.0 * math.sqrt(48))

    def test_longer_solint_raises_the_snr(self):
        assert solution_snr(0.01, 1.0, 27, 60) > solution_snr(0.01, 1.0, 27, 1)

    @pytest.mark.parametrize("arguments", [(1.0, 0.0, 10, 1), (1.0, 1.0, 1, 1), (1.0, 1.0, 10, 0)])
    def test_invalid_inputs_raise(self, arguments):
        with pytest.raises(ValueError):
            solution_snr(*arguments)
