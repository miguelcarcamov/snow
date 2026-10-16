# **SNOW**

## ca**S**a pytho**N** self-calibrati**O**n frame**W**ork

**SNOW** is an object-oriented Python framework for radio-interferometric self-calibration that simplifies the process of writing self-calibration scripts. Instead of rewriting calibration workflows for each dataset, SNOW provides a unified interface to run imagers (CASA tclean, or your own) and self-calibration algorithms (phase, amplitude, amplitude-phase) in a single, streamlined workflow.

### Features

- **Object-oriented design**: Clean, modular architecture for easy customization
- **Multiple imagers supported**: CASA tclean, or your own; Pyralysis (GPU) and WSClean are planned
- **Flexible self-calibration**: Phase-only, amplitude-only, or combined amplitude-phase calibration
- **CASA integration**: Built on CASA (Common Astronomy Software Applications) for radio astronomy data processing
- **Reproducible workflows**: Version-controlled calibration pipelines
- **Docker support**: Pre-built container images for easy deployment

## Requirements

- **Python**: `>=3.12, <3.14`, the versions with both CASA wheels and Pyralysis support
- **Dependencies**: casatools, casatasks (CASA >= 6.7), numpy, astropy and reproject. The
  package declares lower bounds; [`requirements.txt`](requirements.txt) pins the set tested in CI
  and used in the Docker images.
- **CASA data**: CASA downloads its measures data the first time `casatools` is imported
  (see [External Data](https://casadocs.readthedocs.io/en/stable/notebooks/external-data.html)).
  Without network access, point `measurespath` to an existing copy.

On Debian, Ubuntu or Arch, the `manylinux_2_34` casatools wheel that pip picks fails to import
(`OPENSSL_3.0.1 not found`). Install from `requirements.txt`, which uses the working
`manylinux_2_28` wheel, or install that wheel first.

## Installation

### Option 1: From PyPI (Recommended)

```bash
pip install snow-selfcal
```

The package is imported as `snow`. (`pip install snow` installs an unrelated project.) Until
the first release is published, install from GitHub instead.

### Option 2: From GitHub

```bash
pip install -U git+https://github.com/miguelcarcamov/snow.git
```

### Option 3: Using Conda

Create a conda environment from the provided environment file:

```bash
conda env create -f environment.yml
conda activate snow-env
```

The environment provides Python 3.12 and the pinned dependencies from `requirements.txt`; then
install snow itself with `pip install -e .` (from a clone) or `pip install snow-selfcal`.

### Option 4: From Source

```bash
git clone https://github.com/miguelcarcamov/snow
cd snow
pip install .
```

### Option 5: Development Installation

For development work with editable installation:

```bash
git clone https://github.com/miguelcarcamov/snow
cd snow
pip install -e .
```

### Option 6: Docker Container

Pull the pre-built Docker image:

```bash
docker pull ghcr.io/miguelcarcamov/snow:latest
```

The image has Python 3.12, snow with its pinned dependencies, and CASA's measures data and
observatory table in `/opt/casa/data`, readable by any user (also under Apptainer/Singularity). To stay small it has
no CUDA, and it skips the rest of casarundata (antenna responses and simulator files, used by
e.g. the `awproject` gridder); mount a full copy at `/opt/casa/data` if you need it:

```bash
docker run --rm -v "$PWD":/work -w /work ghcr.io/miguelcarcamov/snow:latest python main_ALMA.py target.ms out
```

CASA needs a writable home directory. Apptainer/Singularity provides yours; with `docker run -u`,
also pass `-e HOME=/tmp`.

`snow.utils.reproject()` needs the optional `reproject` package: `pip install "snow-selfcal[reproject]"`.

## Quick Start

### Basic Usage Example

Phase, then amplitude-phase self-calibration of an ALMA dataset. This is
[`main_files/main_ALMA.py`](main_files/main_ALMA.py); a test keeps the two identical.

```python
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
```

Run the script:

```bash
python main_files/main_ALMA.py target.ms target_selfcal
```

How the steps fit together:

- Each self-calibration step continues from the previous one with `previous_selfcal=...`: it
  uses that step's measurement set, model and accepted calibration tables, which stay applied.
  Create it after the previous step's `run()`.
- Amplitude steps (`Ampcal`, `AmpPhasecal`) need phase solutions first; an order that would
  discard or skip them raises a `ValueError` before any data is touched.
- `spwmap` is derived per calibration table from `combine`, so it rarely needs to be set.

More examples for VLA and e-MERLIN are in [`main_files/`](main_files/).

## API Overview

### Imaging Classes

SNOW supports multiple imaging backends:

- **`Tclean`**: CASA tclean task wrapper for standard imaging
- **`WSClean`**: placeholder, not implemented yet
- **`Imager`**: Base class for custom imagers

### Self-Calibration Classes

- **`Phasecal`**: Phase-only self-calibration
- **`Ampcal`**: Amplitude-only self-calibration
- **`AmpPhasecal`**: Combined amplitude-phase self-calibration
- **`Selfcal`**: Base class for self-calibration algorithms

### Utility Modules

- **`snow.utils.image_utils`**: Image processing utilities
- **`snow.utils.selfcal_utils`**: Self-calibration helper functions

## Project Structure

```
snow/
├── src/snow/
│   ├── imaging/          # Imaging backends (tclean; wsclean planned)
│   ├── selfcalibration/  # Self-calibration algorithms
│   └── utils/            # Utility functions
├── main_files/           # Example scripts for different telescopes
├── requirements.txt      # Python dependencies
├── environment.yml       # Conda environment specification
└── pyproject.toml       # Project metadata and build configuration
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the terms specified in the `LICENSE` file.

## Citation

If you use SNOW in your research, please cite:

```bibtex
@software{snow,
  author = {Cárcamo, Miguel},
  title = {SNOW: CASA Python Self-calibration Framework},
  url = {https://github.com/miguelcarcamov/snow},
  version = {<version>},
  year = {<year>}
}
```

## Support

- **Issues**: Report bugs or request features on [GitHub Issues](https://github.com/miguelcarcamov/snow/issues)
- **Documentation**: See the [CASA documentation](https://casadocs.readthedocs.io/) for CASA-specific questions
- **Contact**: <miguel.carcamo@usach.cl>

## Acknowledgments

SNOW is built on top of [CASA](https://casa.nrao.edu/) (Common Astronomy Software Applications), developed by the National Radio Astronomy Observatory (NRAO).
