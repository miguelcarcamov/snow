import sys
import os
import numpy as np
from selfcalframework.selfcal import *
from selfcalframework.imager import *
from flagdata import flagdata

if __name__ == '__main__':
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
    solint_phs = ['32s', '16s']
    #solint_amp = ['1h']
    solint_ap = ['32s']

    phscal = Phasecal(minsnr=3.0, solint=solint_phs, combine="spw", Imager=clean_imager_phs, **shared_vars_dict)

    phs_caltable = phscal.run()

    # ampcal = Ampcal(minsnr=2.0, solint=solint_amp, combine="scan",
    #                selfcal_object=parent_selfcal, input_caltable=phs_caltable)

    #amp_caltable = ampcal.run()

    apcal = AmpPhasecal(minsnr=3.0, solint=solint_ap, combine="", input_caltable=phs_caltable, Imager=clean_imager_ampphs, **shared_vars_dict)

    apcal.run()

    apcal.selfcal_output(overwrite=True)
