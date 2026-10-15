"""
Phase, then amplitude-phase self-calibration of an ALMA dataset imaged with tclean.

Usage: python main_ALMA.py <visfile> <output_prefix> [--plot]
"""
import argparse

from snow.imaging import Tclean
from snow.selfcalibration import AmpPhasecal, Phasecal


def parse_arguments():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("visfile", help="Measurement set with the calibrated data in DATA")
    parser.add_argument("output", help="Prefix of the output images and calibration tables")
    parser.add_argument("--plot", action="store_true", help="Plot the gain solutions")
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_arguments()

    # Automasking values for long baselines: https://casaguides.nrao.edu/index.php/Automasking_Guide
    imaging_parameters = dict(
        inputvis=args.visfile,
        output=args.output,
        niter=100,
        M=1024,
        N=1024,
        cell="0.005arcsec",
        stokes="I",
        data_column="corrected",
        robust=0.5,
        specmode="mfs",
        deconvolver="hogbom",
        gridder="standard",
        save_model=True,
        use_mask="auto-multithresh",
        sidelobe_threshold=3.0,
        noise_threshold=5.0,
        min_beam_frac=0.3,
        low_noise_threshold=1.5,
        negative_threshold=0.0,
        interactive=False  # True to edit the masks by hand
    )
    clean_imager_phs = Tclean(threshold="0.1mJy", **imaging_parameters)
    clean_imager_ampphs = Tclean(threshold="0.025mJy", **imaging_parameters)

    # Parameters shared by every self-calibration step. spwmap is derived from combine.
    shared_parameters = dict(
        minsnr=3.0, minblperant=6, refant="DA51", gaintype="T", want_plot=args.plot
    )

    phscal = Phasecal(
        visfile=args.visfile,
        imager=clean_imager_phs,
        solint=["inf", "600s"],
        combine="spw",
        **shared_parameters
    )
    phscal.run()

    # Continues from the phase self-calibration (its measurement set, model and tables), so it is
    # created after phscal.run()
    apcal = AmpPhasecal(
        previous_selfcal=phscal,
        imager=clean_imager_ampphs,
        solint=["inf"],
        combine="",
        **shared_parameters
    )
    apcal.run()

    apcal.selfcal_output(overwrite=True)
