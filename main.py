from selfcal import *
from imager import *
import sys
import numpy as np
sys.path.append('./')

if __name__ == '__main__':
    visfile = sys.argv[3]
    output = sys.argv[4]
    want_plot = eval(sys.argv[5])

    imager_obj = Imager(inputvis=visfile, output=output,
                        niter=100, M=1024, N=1024, deltax="0.02arcsec", stokes="I", datacolumn="corrected", robust=0.5)

    clean_imager = Clean(specmode="mfs", deconvolver="hogbom", gridder="standard",
                         pbcor=True, savemodel="modelcolumn", imager_object=imager_obj, interactive=True)

    parent_selfcal = Selfcal(visfile=clean_imager.inputvis,
                             imagename=clean_imager.output, minblperant=2, refant="Kn,Mk2,Cm,Pi,De,Da", spwmap=[0, 0, 0, 0, 0, 0, 0, 0], Imager=clean_imager, want_plot=want_plot)

    solint_phs = ['128s', '64s', '32s', '16s']
    solint_amp = ['1h']
    solint_ap = ['inf']

    phscal = Phasecal(minsnr=2.0, solint=solint_phs,
                      combine="spw", selfcal_object=parent_selfcal)

    phs_caltable = phscal.run()

    ampcal = Ampcal(minsnr=2.0, solint=solint_amp, combine="scan",
                    selfcal_object=parent_selfcal, input_caltable=phs_caltable)

    amp_caltable = ampcal.run()

    apcal = AmpPhasecal(minsnr=2.0, solint=solint_ap, combine="",
                        selfcal_object=parent_selfcal, input_caltable=amp_caltable)

    apcal.run()

    parent_selfcal.selfcal_output()
