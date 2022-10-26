# **SNOW**

## ca**S**a pytho**N** self-calibrati**O**n frame**W**ork

Many radio-astronomers repeat the process of writing different scripts for self-calibration
depending on their datasets. This repository holds an object-oriented Framework for self-calibration
of radio-interferometric datasets that will help radio astronomers to minimize the tedious work of
writing self-calibration scripts once again. The idea is to call just one main Python script that
will run an imager (tclean, wsclean, gpuvmem, rascil, etc.) and one or multiple self-calibration
objects (phase, amplitude, amplitude-phase) having the self-calibrated dataset as a result.

## Requirements

1. `Python == 3.8`
2. Check CASA pip current version requirements [here](https://casadocs.readthedocs.io/en/stable/notebooks/introduction.html#Modular-Packages).
3. Check the `requirements.txt` file.

## Installation

### From PYPI repository

- `pip install snow`

### From Github

- `pip install -U git+https://github.com/miguelcarcamov/snow`

### From source

```bash
git clone https://github.com/miguelcarcamov/snow
cd snow
pip install .
```

### From source as developer

```bash
git clone https://github.com/miguelcarcamov/snow
cd snow
pip install e .
```

## Using docker container

```bash
docker pull ghcr.io/miguelcarcamov/snow:latest
```

## Run snow

```Python
# Import the modules that you want to use
import sys
from snow.selfcalibration import Phasecal, AmpPhasecal
from snow.imaging import Tclean

if __name__ == '__main__':
 # This step is up to you, and option to capture your arguments from terminal is using sys.argv
 visfile = sys.argv[3]
 output = sys.argv[4]
 want_plot = eval(sys.argv[5])

 # Table for automasking on long or short baselines can be found here: https://casaguides.nrao.edu/index.php/Automasking_Guide
 # The default clean object will use automasking values for short baselines
 # In this case we will use automasking values for long baselines
 # Create different imagers with different thresholds (this is optional, you can create just one)
 clean_imager_phs = Tclean(inputvis=visfile, output=output, niter=100, M=1024, N=1024, cell="0.005arcsec",
                           stokes="I", datacolumn="corrected", robust=0.5,
                           specmode="mfs", deconvolver="hogbom", gridder="standard",
                           savemodel=True, usemask='auto-multithresh', threshold="0.1mJy", sidelobethreshold=3.0,
                           noisethreshold=5.0,
                           minbeamfrac=0.3, lownoisethreshold=1.5, negativethreshold=0.0, interactive=True)

 clean_imager_ampphs = Tclean(inputvis=visfile, output=output, niter=100, M=1024, N=1024, cell="0.005arcsec",
                              stokes="I", datacolumn="corrected", robust=0.5,
                              specmode="mfs", deconvolver="hogbom", gridder="standard",
                              savemodel=True, usemask='auto-multithresh', threshold="0.025mJy",
                              sidelobethreshold=3.0,
                              noisethreshold=5.0,
                              minbeamfrac=0.3, lownoisethreshold=1.5, negativethreshold=0.0, interactive=True)

 # This is a dictionary with shared variables between self-cal objects
 shared_vars_dict = {'visfile': visfile, 'minblperant': 6, 'refant': "DA51", 'spwmap': [
  0, 0, 0, 0], 'gaintype': 'T', 'want_plot': want_plot}

 # Create your solution intervals
 solint_phs = ['inf', '600s']
 solint_ap = ['inf']

 # Create your phasecal object
 phscal = Phasecal(minsnr=3.0, solint=solint_phs, combine="spw", imager=clean_imager_phs, **shared_vars_dict)
 # Run it!
 phscal.run()

 # If we are happy with the result of the only-phase self-cal we can end the code here, if not...
 # Create the amplitude-phase self-cal object
 apcal = AmpPhasecal(minsnr=3.0, solint=solint_ap, combine="", previous_selfcal=phscal, imager=clean_imager_ampphs,
                     **shared_vars_dict)
 # Run it
 apcal.run()
 # Get your splitted final MS
 apcal.selfcal_output(overwrite=True)
```

Then you can simply run the main script using `python yourscript.py <arguments>`
