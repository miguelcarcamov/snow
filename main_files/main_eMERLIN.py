from selfcalframework.selfcal import *
from selfcalframework.imager import *
from casatasks import flagdata
from casatasks import mstransform
import sys
import numpy as np
import os

if __name__ == '__main__':

    print(sys.argv)
    visfile = sys.argv[1]
    output = sys.argv[2]
    want_plot = eval(sys.argv[3])

    selfcal_vis = visfile.split("./../")[1][:-3]+".flagged.ms"
    if not os.path.exists(selfcal_vis): mstransform(vis=visfile, outputvis=selfcal_vis, datacolumn="corrected")

    clean_imager = Clean(inputvis=selfcal_vis, output=output, niter=10000, M=1024, N=1024, cell="0.02arcsec",
                         stokes="I", datacolumn="corrected", robust=2.0, scales=[0, 3, 5, 10 ,15, 20, 30, 40, 50, 80], specmode="mfs", deconvolver="multiscale", gridder="standard",
                         pbcor=True, savemodel=True, nsigma=4.0, interactive=False, cycleniter=100, usemask='auto-multithresh', sidelobethreshold=1.0, noisethreshold=8.0,
                         minbeamfrac=0.2, lownoisethreshold=1.5, negativethreshold=0.0)

    shared_vars_dict = {'visfile': clean_imager.getVis(), 'minblperant': 2, 'refant': "Kn,Cm,Mk2,Pi,De,Da", 'spwmap': [0, 0, 0, 0, 0, 0, 0, 0], 'gaintype': 'G', 'want_plot': want_plot}

    #solint_phs = ['128s', '64s', '32s', '16s']
    solint_phs = ['inf', '3min', '2min', '1min', '30s', '15s']
    varchange_phs = {'nsigma' : [3.0, 2.0, 2.0, 1.5, 1.5, 1.5]}
    solint_amp = ['1h']
    solint_ap = ['inf']

    phscal = Phasecal(minsnr=2.0, solint=solint_phs,
                      combine="spw", varchange=varchange_phs, Imager=clean_imager, **shared_vars_dict)

    # Backup MS to the state before self-cal
    flagmanager(vis=selfcal_vis, mode='save', versionname='before_selfcal_flagging')

    # Flag residual RFI
    mode="rflag"
    flagdata(vis=selfcal_vis, datacolumn="data", spw="0", correlation="LL,RR", mode=mode, action="apply", flagbackup=False)
    flagdata(vis=selfcal_vis, datacolumn="data", spw="0", mode="extend", extendflags=False, action="apply", flagbackup=False)
    flagdata(vis=selfcal_vis, datacolumn="data", spw="1~2", correlation="LL,RR", mode=mode, action="apply", flagbackup=False)
    flagdata(vis=selfcal_vis, datacolumn="data", spw="1~2", mode="extend", extendflags=False, action="apply", flagbackup=False)
    flagdata(vis=selfcal_vis, datacolumn="data", spw="3~5", correlation="LL,RR", mode=mode, action="apply", flagbackup=False)
    flagdata(vis=selfcal_vis, datacolumn="data", spw="3~5", mode="extend", extendflags=False, action="apply", flagbackup=False)
    flagdata(vis=selfcal_vis, datacolumn="data", spw="5~7", correlation="LL,RR", mode=mode, action="apply", flagbackup=False)
    flagdata(vis=selfcal_vis, datacolumn="data", spw="5~7", mode="extend", extendflags=False, action="apply", flagbackup=False)

    # Backup MS to the state after self-cal
    flagmanager(vis=selfcal_vis, mode='save', versionname='after_selfcal_flagging')

    # Run phase-loop self-calibration
    phscal.run()

    # ampcal = Ampcal(minsnr=2.0, solint=solint_amp, combine="scan",
    #                selfcal_object=parent_selfcal, input_caltable=phs_caltable)

    #amp_caltable = ampcal.run()

    #apcal = AmpPhasecal(minsnr=2.0,
    #                    solint=solint_ap, combine="", input_caltable=phs_caltable, Imager=clean_imager, **shared_vars_dict)

    #apcal.run()

    phscal.selfcal_output(overwrite=True)
