"""
Amplitude-phase self-calibration of an e-MERLIN target, on top of existing phase solutions, imaged
with multiscale tclean. Optionally flags residual RFI before self-calibrating.

Usage: python main_eMERLIN.py <visfile> <output_prefix> <field> [--phase-caltable pcal1] [--flag]
       [--plot]
"""
import argparse
import os
from pathlib import Path

from casatasks import flagdata, flagmanager, mstransform

from snow.imaging import Tclean
from snow.selfcalibration import AmpPhasecal

# Spectral window groups flagged together with rflag
RFLAG_SPW_GROUPS = ("0", "1~2", "3~5", "5~7")


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("visfile", help="Calibrated measurement set (CORRECTED_DATA is used)")
    parser.add_argument("output", help="Prefix of the output images and calibration tables")
    parser.add_argument("field", help="Target field to self-calibrate")
    parser.add_argument(
        "--phase-caltable",
        default="pcal1",
        help="Existing phase calibration table to pre-apply (solutions combined over spws)"
    )
    parser.add_argument(
        "--flag", action="store_true", help="Flag residual RFI with rflag before self-calibrating"
    )
    parser.add_argument("--plot", action="store_true", help="Plot the gain solutions")
    return parser.parse_args()


def flag_residual_rfi(vis: str) -> None:
    flagmanager(vis=vis, mode='save', versionname='before_selfcal_flagging')
    for spw in RFLAG_SPW_GROUPS:
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
        flag_residual_rfi(selfcal_vis)

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
        nsigma=3.0,
        interactive=False,
        cycle_niter=100,
        use_mask="user",
        mask=args.output + "_ph2.mask"
    )

    # Imaging parameters that change on each solution interval
    varchange_imager = {
        "nsigma": [2.0, 2.0, 2.0],
        "use_mask": ["user", "user", "auto-multithresh"],
        "mask": [args.output + "_ph2.mask", args.output + "_ph2.mask", ""],
    }

    apcal = AmpPhasecal(
        visfile=selfcal_vis,
        imager=clean_imager,
        input_caltable=args.phase_caltable,
        # Map of the input phase table, whose solutions are combined over the 8 spws
        spwmap=[0, 0, 0, 0, 0, 0, 0, 0],
        solint=["14min", "7min", "3.5min"],
        varchange_imager=varchange_imager,
        combine="spw",
        applymode="calonly",
        uvrange=">60klambda",
        minsnr=3.0,
        minblperant=2,
        refant="Kn,Cm,Mk2,Pi,De,Da",
        gaintype="G",
        flag_dataset=True,
        restore_psnr=True,
        want_plot=args.plot
    )
    apcal.run()

    apcal.selfcal_output(overwrite=True)
