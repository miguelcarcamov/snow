# **SNOW**

## ca**S**a pytho**N** self-calibrati**O**n frame**W**ork

**SNOW** is an object-oriented Python framework for radio-interferometric self-calibration that simplifies the process of writing self-calibration scripts. Instead of rewriting calibration workflows for each dataset, SNOW provides a unified interface to run imagers (tclean, wsclean, gpuvmem, rascil, etc.) and self-calibration algorithms (phase, amplitude, amplitude-phase) in a single, streamlined workflow.

### Features

- **Object-oriented design**: Clean, modular architecture for easy customization
- **Multiple imagers supported**: CASA tclean, wsclean, GPUvmem, and more
- **Flexible self-calibration**: Phase-only, amplitude-only, or combined amplitude-phase calibration
- **CASA integration**: Built on CASA (Common Astronomy Software Applications) for radio astronomy data processing
- **Reproducible workflows**: Version-controlled calibration pipelines
- **Docker support**: Pre-built container images for easy deployment

## Requirements

- **Python**: `>=3.10, <3.11` (Python 3.10 only, due to CASA compatibility requirements)
- **CASA packages**: See [CASA pip package requirements](https://casadocs.readthedocs.io/en/stable/notebooks/introduction.html#Modular-Packages) for current version requirements
- **Dependencies**: See `requirements.txt` for a complete list of required packages

### Key Dependencies

- **CASA packages**: casatasks, casatools, casashell, and related CASA modules
- **Scientific Python**: numpy (1.26.0), scipy (1.15.0), astropy (6.1.0)
- **Astronomy tools**: reproject (0.14.1), python-casacore (3.7.1)

## Installation

### Option 1: From PyPI (Recommended)

```bash
pip install snow
```

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

The environment file includes all system dependencies and will automatically install Python packages from `requirements.txt`.

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

The Docker image includes all dependencies and is ready to use. See the [Docker documentation](https://docs.docker.com/) for usage instructions.

## Quick Start

### Basic Usage Example

```python
import sys
from snow.selfcalibration import Phasecal, AmpPhasecal
from snow.imaging import Tclean

if __name__ == '__main__':
    # Get input parameters (adjust as needed for your use case)
    visfile = sys.argv[1]  # Input measurement set
    output = sys.argv[2]   # Output prefix
    want_plot = eval(sys.argv[3])  # Whether to generate plots

    # Create imager for phase calibration
    # See CASA automasking guide: https://casaguides.nrao.edu/index.php/Automasking_Guide
    clean_imager_phs = Tclean(
        inputvis=visfile,
        output=output,
        niter=100,
        M=1024,
        N=1024,
        cell="0.005arcsec",
        stokes="I",
        datacolumn="corrected",
        robust=0.5,
        specmode="mfs",
        deconvolver="hogbom",
        gridder="standard",
        savemodel=True,
        usemask='auto-multithresh',
        threshold="0.1mJy",
        sidelobethreshold=3.0,
        noisethreshold=5.0,
        minbeamfrac=0.3,
        lownoisethreshold=1.5,
        negativethreshold=0.0,
        interactive=True
    )

    # Create imager for amplitude-phase calibration (with lower threshold)
    clean_imager_ampphs = Tclean(
        inputvis=visfile,
        output=output,
        niter=100,
        M=1024,
        N=1024,
        cell="0.005arcsec",
        stokes="I",
        datacolumn="corrected",
        robust=0.5,
        specmode="mfs",
        deconvolver="hogbom",
        gridder="standard",
        savemodel=True,
        usemask='auto-multithresh',
        threshold="0.025mJy",
        sidelobethreshold=3.0,
        noisethreshold=5.0,
        minbeamfrac=0.3,
        lownoisethreshold=1.5,
        negativethreshold=0.0,
        interactive=True
    )

    # Shared parameters for self-calibration objects
    shared_vars_dict = {
        'visfile': visfile,
        'minblperant': 6,
        'refant': "DA51",
        'spwmap': [0, 0, 0, 0],
        'gaintype': 'T',
        'want_plot': want_plot
    }

    # Solution intervals for phase calibration
    solint_phs = ['inf', '600s']

    # Solution intervals for amplitude-phase calibration
    solint_ap = ['inf']

    # Phase-only self-calibration
    phscal = Phasecal(
        minsnr=3.0,
        solint=solint_phs,
        combine="spw",
        imager=clean_imager_phs,
        **shared_vars_dict
    )
    phscal.run()

    # Amplitude-phase self-calibration (optional, only if phase calibration is successful)
    apcal = AmpPhasecal(
        minsnr=3.0,
        solint=solint_ap,
        combine="",
        previous_selfcal=phscal,
        imager=clean_imager_ampphs,
        **shared_vars_dict
    )
    apcal.run()

    # Output the final calibrated measurement set
    apcal.selfcal_output(overwrite=True)
```

Run the script:

```bash
python yourscript.py <visfile> <output_prefix> True
```

## API Overview

### Imaging Classes

SNOW supports multiple imaging backends:

- **`Tclean`**: CASA tclean task wrapper for standard imaging
- **`WSClean`**: WSClean imager interface
- **`GPUvmem`**: GPU-accelerated imaging with GPUvmem
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
│   ├── imaging/          # Imaging backends (tclean, wsclean, gpuvmem)
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
