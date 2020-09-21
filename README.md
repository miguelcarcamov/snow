# Object Oriented Framework for Self-calibration of radio-interferometric datasets

Many radioastronomers repeat the process of writing different scripts for self-calibration depending on their datasets. This repository holds an object oriented Framework for self-calibration of radio-interferometric datasets that will help radioastronomers to minimize the tedious work of writing self-calibration scripts once again. The idea is to call just one main Python script that will run an imager (tclean, wsclean, gpuvmem, rascil, etc.) and one or multuple self-calibration objects (phase, amplitude, amplitude-phase) having the self-calibrated dataset as a result. **It is important to recall that this repository is heavily under development!**
## Requirements

1. CASA (https://casa.nrao.edu/casa_obtaining.shtml)

## Installation
We need to install the selfcalframework modules in CASA in order to call the different objects (selfcal and imager). The installation is very similar to the astropy installation in CASA.

- If you want to modify or develop modules and test them:

  1. Open CASA in the repository folder of the framework
  2. Install pip inside CASA
  ```Python
  CASA <2>: from setuptools.command import easy_install
  CASA <3>: easy_install.main(['--user', 'pip'])
  ```
  3. Quit CASA, re-open it and install the selfcalframework modules
  ```Python
  CASA <2>: import subprocess, sys
  CASA <3>: subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--user', '-e', '.'])
  ```
  4. Then close CASA again and open it, and you should be able to import `selfcalframework` from CASA or your CASA scripts
  ```Python
  CASA <2>: from selfcalframework.imager import *
  CASA <3>: from selfcalframework.selfcal import *
  ```
- If you just want to use the modules inside CASA:

  1. Open CASA
  2. Install pip inside CASA
  ```Python
  CASA <2>: from setuptools.command import easy_install
  CASA <3>: easy_install.main(['--user', 'pip'])
  ```
  3. Quit CASA, re-open it and install the selfcalframework modules
  ```Python
  CASA <2>: import subprocess, sys
  CASA <3>: subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--user', 'selfcalframework'])
  ```
  4. Then close CASA again and open it, and you should be able to import `selfcalframework` from CASA or your CASA scripts
  ```Python
  CASA <2>: from selfcalframework.imager import *
  CASA <3>: from selfcalframework.selfcal import *
  ```


## Run your scripts

In the `main_files` folder there is a set of examples script to run your self-calibration. As a example we will show one of them here:

```Python
# Import the modules that you want to use
import sys
import numpy as np
from selfcalframework.selfcal import *
from selfcalframework.imager import *

if __name__ == '__main__':
    # This step is up to you, and option to capture your arguments from terminal is using sys.argv
    visfile = sys.argv[3]
    output = sys.argv[4]
    want_plot = eval(sys.argv[5])

    # Table for automasking on long or short baselines can be found here: https://casaguides.nrao.edu/index.php/Automasking_Guide
      # The default clean object will use automasking values for short baselines
      # In this case we will use automasking values for long baselines
      clean_imager_phs = Clean(inputvis=visfile, output=output, niter=100, M=1024, N=1024, cell="0.005arcsec",
                               stokes="I", datacolumn="corrected", robust=0.5,
                               specmode="mfs", deconvolver="hogbom", gridder="standard",
                               savemodel=True, usemask='auto-multithresh', threshold="0.1mJy", sidelobethreshold=3.0, noisethreshold=5.0,
                               minbeamfrac=0.3, lownoisethreshold=1.5, negativethreshold=0.0, interactive=True)

      clean_imager_ampphs = Clean(inputvis=visfile, output=output, niter=100, M=1024, N=1024, cell="0.005arcsec",
                                  stokes="I", datacolumn="corrected", robust=0.5,
                                  specmode="mfs", deconvolver="hogbom", gridder="standard",
                                  savemodel=True, usemask='auto-multithresh', threshold="0.025mJy", sidelobethreshold=3.0, noisethreshold=5.0,
                                  minbeamfrac=0.3, lownoisethreshold=1.5, negativethreshold=0.0, interactive=True)

      shared_vars_dict = {'visfile': visfile, 'minblperant': 6, 'refant': "DA51", 'spwmap': [
          0, 0, 0, 0], 'gaintype': 'T', 'want_plot': want_plot}

      #solint_phs = ['128s', '64s', '32s', '16s']
      solint_phs = ['600s', 'inf']
      #solint_amp = ['1h']
      solint_ap = ['inf']

      phscal = Phasecal(minsnr=3.0, solint=solint_phs, combine="spw", Imager=clean_imager_phs, **shared_vars_dict)

      phscal.run()

      # ampcal = Ampcal(minsnr=2.0, solint=solint_amp, combine="scan",
      #                selfcal_object=parent_selfcal, input_caltable=phs_caltable)

      #amp_caltable = ampcal.run()

      apcal = AmpPhasecal(minsnr=3.0, solint=solint_ap, combine="", selfcal_object=phscal, Imager=clean_imager_ampphs, **shared_vars_dict)

      apcal.run()

      apcal.selfcal_output(overwrite=True)
```

Then you can simply run the main script using `casa -c yourscript.py <arguments>`
