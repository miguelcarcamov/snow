from selfcalframework.selfcal import *
from selfcalframework.imager import *
import sys
import numpy as np

if __name__ == '__main__':
    visfile = sys.argv[3]
    output = sys.argv[4]
    want_plot = eval(sys.argv[5])

    clean_imager = Clean(inputvis=visfile, output=output, niter=100, M=1024, N=1024, cell="0.02arcsec",
                         stokes="I", datacolumn="corrected", robust=0.5, specmode="mfs", deconvolver="hogbom", gridder="standard",
                         pbcor=True, savemodel=True, interactive=True)

    shared_vars_dict = {'visfile': clean_imager.getVis(), 'minblperant': 2, 'refant': "Kn,Mk2,Cm,Pi,De,Da", 'spwmap': [0, 0, 0, 0, 0, 0, 0, 0], 'gaintype': 'G', 'want_plot': want_plot}

    solint_phs = ['128s', '64s', '32s', '16s']
    solint_amp = ['1h']
    solint_ap = ['inf']

    phscal = Phasecal(minsnr=2.0, solint=solint_phs,
                      combine="spw", Imager=clean_imager, **shared_vars_dict)

    phscal.run()

    # ampcal = Ampcal(minsnr=2.0, solint=solint_amp, combine="scan",
    #                selfcal_object=parent_selfcal, input_caltable=phs_caltable)

    #amp_caltable = ampcal.run()

    #apcal = AmpPhasecal(minsnr=2.0,
    #                    solint=solint_ap, combine="", input_caltable=phs_caltable, Imager=clean_imager, **shared_vars_dict)

    #apcal.run()

    phscal.selfcal_output(overwrite=True)
