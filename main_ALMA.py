import sys
import os
import numpy as np
sys.path.append('/home/miguel.carcamo/SELFCAL_TEST/objectoriented_selfcal')
print os.path.abspath('.')
from selfcal import *
from imager import *

if __name__ == '__main__':
    visfile = sys.argv[3]
    output = sys.argv[4]
    want_plot = eval(sys.argv[5])

    imager_obj = Imager(inputvis=visfile, output=output,
                        niter=100, M=1024, N=1024, deltax="0.005arcsec", stokes="I", datacolumn="corrected", robust=0.5)

    # Table for long or short baselines can be found here: https://casaguides.nrao.edu/index.php/Automasking_Guide
    # The default clean object will use automasking values for short baselines
    # In this case we will use automasking values for long baselines
    clean_imager = Clean(specmode="mfs", deconvolver="hogbom", gridder="standard",
                         pbcor=True, savemodel="modelcolumn", imager_object=imager_obj,
                         usemask='auto-multithresh', sidelobethreshold=3.0, noisethreshold=5.0,
                         minbeamfrac=0.3, lownoisethreshold=1.5, negativethreshold=0.0, interactive=True)

    parent_selfcal = Selfcal(visfile=clean_imager.inputvis,
                             imagename=clean_imager.output, minblperant=6, refant="DV03", spwmap=[0, 0, 0, 0], Imager=clean_imager, want_plot=want_plot)

    #solint_phs = ['128s', '64s', '32s', '16s']
    solint_phs = ['32s', '16s']
    #solint_amp = ['1h']
    solint_ap = ['32s']

    phscal = Phasecal(minsnr=2.0, solint=solint_phs,
                      combine="spw", selfcal_object=parent_selfcal)

    phs_caltable = phscal.run()

    #ampcal = Ampcal(minsnr=2.0, solint=solint_amp, combine="scan",
    #                selfcal_object=parent_selfcal, input_caltable=phs_caltable)

    #amp_caltable = ampcal.run()

    #apcal = AmpPhasecal(minsnr=2.0, solint=solint_ap, combine="",
    #                    selfcal_object=parent_selfcal, input_caltable=phs_caltable)

    #apcal.run()

    parent_selfcal.selfcal_output(overwrite=True)
