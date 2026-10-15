"""
Phase, then amplitude-phase self-calibration of a VLA dataset imaged with tclean.

Usage: python main_VLA.py <visfile> <output_prefix> [--plot]
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
    clean_imager = Tclean(
        inputvis=args.visfile,
        output=args.output,
        niter=100,
        M=1024,
        N=1024,
        cell="0.3arcsec",
        stokes="I",
        data_column="corrected",
        robust=0.5,
        specmode="mfs",
        deconvolver="hogbom",
        gridder="standard",
        pbcor=True,
        save_model=True,
        use_mask="auto-multithresh",
        sidelobe_threshold=1.25,
        noise_threshold=5.0,
        min_beam_frac=0.1,
        low_noise_threshold=2.0,
        negative_threshold=0.0,
        interactive=False  # True to edit the masks by hand
    )

    # Parameters shared by every self-calibration step. spwmap is derived from combine.
    shared_parameters = dict(
        minsnr=3.0, minblperant=4, refant="VA05", gaintype="G", want_plot=args.plot
    )

    phscal = Phasecal(
        visfile=args.visfile,
        imager=clean_imager,
        solint=["inf"],
        combine="spw",
        **shared_parameters
    )
    phscal.run()

    # Continues from the phase self-calibration (its measurement set, model and tables), so it is
    # created after phscal.run()
    apcal = AmpPhasecal(
        previous_selfcal=phscal,
        imager=clean_imager,
        solint=["inf"],
        combine="",
        **shared_parameters
    )
    apcal.run()

    apcal.selfcal_output(overwrite=True)
