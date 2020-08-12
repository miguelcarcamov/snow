from selfcalframework.selfcal import *
from selfcalframework.imager import *
import sys
import numpy as np
from concat import concat
from initweights import initweights
from statwt import statwt
from selfcalframework.selfcal_utils import *

if __name__ == '__main__':
    visfiles = sys.argv[3].split(',')
    print("Visfiles: ", visfiles)
    output = sys.argv[4]
    want_plot = eval(sys.argv[5])
    normweights = eval(sys.argv[6])

    imager_objects = []
    clean_objects = []
    selfcal_objects = []
    outputs = [output+"/D_img", output+"/DC_img", output+"/DCBA_img"]
    configs = ["D", "C", "BA"]
    visnames = []
    visnames.append(visfiles[0])
    concatvis = []
    concatvis.append(visfiles[0])
    concatname = "D"
    #solint_phs = ['128s', '64s', '32s', '16s']
    #solint_phs = ['32s', '16s']
    solint_phs = ['inf']
    #solint_amp = ['1h']
    solint_ap = ['inf']

    deltax_vector = ["0.3arcsec", "0.3arcsec", "0.4arcsec"]

    if normweights:
        for i in visfiles:
            print("Normalizing weights in dataset:", i)
            initweights(vis=i, wtmode='nyq')

    for i in range(0, len(visfiles)):

        if(i > 0):
            concatname += "+" +configs[i]
            visnames.append(concatname + ".ms")
            concatvis.append(visfiles[i])
            print("Selfcal is going to concat: ", concatvis)
            concat(concatvis=visnames[i], vis=concatvis)
            concatvis=[]
            concatvis.append(visnames[i])
        print("Concat name:", concatname)
        print("Visnames: ", visnames)
        spwmap = [0] * getTableRows(visnames[i] + '/SPECTRAL_WINDOW')

        clean_imager = Clean(inputvis=visnames[i], output=outputs[i], niter=100, M=1024, N=1024, cell=deltax_vector[i], stokes="I", datacolumn="corrected",
                                 robust=0.5, specmode="mfs", deconvolver="hogbom", gridder="standard",
                                 pbcor=True, savemodel="modelcolumn", usemask='auto-multithresh', sidelobethreshold=1.25, noisethreshold=5.0,
                                 minbeamfrac=0.1, lownoisethreshold=2.0, negativethreshold=0.0, interactive=True)

        parent_selfcal = Selfcal(visfile=clean_imager.getVis(), minblperant=4, refant="VA05", spwmap=spwmap, Imager=clean_imager, want_plot=want_plot)

        phscal = Phasecal(minsnr=2.0, solint=solint_phs,
                          combine="spw", selfcal_object=parent_selfcal)

        phs_caltable = phscal.run()

        #ampcal=Ampcal(minsnr=2.0, solint=solint_amp, combine="scan",
        #                selfcal_object=parent_selfcal, input_caltable=phs_caltable)

        #amp_caltable=ampcal.run()

        apcal = AmpPhasecal(minsnr=2.0, solint=solint_ap, combine="",
                            selfcal_object=parent_selfcal, input_caltable=phs_caltable)

        apcal.run()

        if(i == len(visfiles)-1):
            parent_selfcal.selfcal_output(overwrite=True)
