"""
Amplitude-phase self-calibration of an e-MERLIN target, on top of existing phase solutions, imaged
with multiscale tclean. Optionally flags residual RFI before self-calibrating.

Usage: python main_eMERLIN.py <visfile> <output_prefix> <field> [--phase-caltable pcal]
       [--phase-caltable-per-spw] [--mask MASK] [--uvrange UVRANGE] [--refant REFANT]
       [--flag [--rflag-spw SPW ...]] [--plot]
"""
import argparse
import os
from pathlib import Path

from casatasks import flagdata, flagmanager, mstransform

from snow.imaging import Tclean
from snow.selfcalibration import AmpPhasecal
from snow.selfcalibration.caltable import derive_spwmap
from snow.utils import count_spectral_windows, get_spw_ids

# Amplitude-phase solution intervals, from long to short
SOLINTS = ["14min", "7min", "3.5min"]


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("visfile", help="Calibrated measurement set (CORRECTED_DATA is used)")
    parser.add_argument("output", help="Prefix of the output images and calibration tables")
    parser.add_argument("field", help="Target field to self-calibrate")
    parser.add_argument(
        "--phase-caltable", default="pcal", help="Existing phase calibration table to pre-apply"
    )
    parser.add_argument(
        "--phase-caltable-per-spw",
        action="store_true",
        help="The phase table has one solution per spw (default: combined over spws)"
    )
    parser.add_argument("--mask", help="CASA mask image; automasking is used if omitted")
    parser.add_argument("--uvrange", default="", help="uv-range used to solve, e.g. '>60klambda'")
    parser.add_argument(
        "--refant", default="Kn,Cm,Mk2,Pi,De,Da", help="Reference antennas, in order of preference"
    )
    parser.add_argument(
        "--flag", action="store_true", help="Flag residual RFI with rflag before self-calibrating"
    )
    parser.add_argument(
        "--rflag-spw",
        action="append",
        default=None,
        help="Spectral window selection flagged together; repeat for several groups "
        "(default: all spws in one pass)"
    )
    parser.add_argument("--plot", action="store_true", help="Plot the gain solutions")
    return parser.parse_args()


def flag_residual_rfi(vis: str, spw_groups: list) -> None:
    flagmanager(vis=vis, mode='save', versionname='before_selfcal_flagging')
    for spw in spw_groups:
        flagdata(
            vis=vis,
            datacolumn="data",
            spw=spw,
            correlation="LL,RR",
            mode="rflag",
            action="apply",
            flagbackup=False
        )
        flagdata(
            vis=vis,
            datacolumn="data",
            spw=spw,
            mode="extend",
            extendflags=False,
            action="apply",
            flagbackup=False
        )
    flagmanager(vis=vis, mode='save', versionname='after_selfcal_flagging')


if __name__ == '__main__':
    args = parse_arguments()

    # Self-calibration solves against DATA, so the calibrated target is split into DATA first
    selfcal_vis = Path(args.visfile).stem + ".flagged.ms"
    if not os.path.exists(selfcal_vis):
        mstransform(
            vis=args.visfile, outputvis=selfcal_vis, datacolumn="corrected", field=args.field
        )

    # Flag before creating the self-calibration object, which works on a copy of selfcal_vis
    if args.flag:
        flag_residual_rfi(selfcal_vis, args.rflag_spw or [""])

    # Cell size and scales suit e-MERLIN at L-band; adjust them to your band and target
    clean_imager = Tclean(
        inputvis=selfcal_vis,
        output=args.output,
        niter=10000,
        M=1024,
        N=1024,
        cell="0.02arcsec",
        stokes="I",
        data_column="corrected",
        robust=2.0,
        scales=[0, 3, 5, 10, 15, 20, 30, 40, 50, 80, 100],
        specmode="mfs",
        deconvolver="multiscale",
        gridder="standard",
        pbcor=False,
        save_model=True,
        interactive=False,
        cycle_niter=100,
        use_mask="user" if args.mask else "auto-multithresh",
        mask=args.mask or ""
    )

    # Clean deeper as the calibration improves: one value per solution interval
    varchange_imager = {"nsigma": [3.0, 2.5, 2.0]}

    # spwmap of the input phase table: solutions combined over spws are stored under the lowest one
    phase_spwmap = derive_spwmap(
        "" if args.phase_caltable_per_spw else "spw", get_spw_ids(selfcal_vis),
        count_spectral_windows(selfcal_vis)
    )

    apcal = AmpPhasecal(
        visfile=selfcal_vis,
        imager=clean_imager,
        input_caltable=args.phase_caltable,
        spwmap=phase_spwmap,
        solint=SOLINTS,
        varchange_imager=varchange_imager,
        combine="spw",
        applymode="calonly",
        uvrange=args.uvrange,
        minsnr=3.0,
        minblperant=2,
        refant=args.refant,
        gaintype="G",
        flag_dataset=True,
        restore_psnr=True,
        want_plot=args.plot
    )
    apcal.run()

    apcal.selfcal_output(overwrite=True)
