from selfcalframework.selfcal import *
from selfcalframework.imager import *
import sys
import numpy as np
from concat import concat
from selfcalframework.selfcal_utils import *

if __name__ == '__main__':
    visfiles = sys.argv[3].split(',')
    output = sys.argv[4]
    want_plot = eval(sys.argv[5])

    first_data = visfiles[0]

    imager_objects = []
    clean_objects = []
    selfcal_objects = []
    outputs = ["D_img", "DC_img", "DCB_img", output]
    arrays = ["B", "C", "D"]
    vis_imaging = []
    vis_imaging.append(first_data)
    visnames = []
    visnames.append(first_data)

    #solint_phs = ['128s', '64s', '32s', '16s']
    solint_phs = ['32s', '16s']
    solint_amp = ['1h']
    solint_ap = ['inf']

    deltax_vector = ["0.4arcsec", "0.3arcsec", "0.3arcsec", "0.3arcsec"]

    for i in range(0, len(visfiles)):

        if(i > 0):
            vis_imaging.append(vis_imaging[i - 1] + arrays[i])
            visnames.append(visfiles[i])
            concat(concatvis=vis_imaging[i], vis=visnames)

        spwmap = [0] * getTableRows(vis_imaging[i] + '/SPECTRAL_WINDOW')

        imager_objects.append(Imager(inputvis=vis_imaging[i], output=outputs[i], niter=100, M=1024,
                                     N=1024, deltax=deltax_vector[i], stokes="I", datacolumn="corrected", robust=0.0))

        clean_objects.append(Clean(specmode="mfs", deconvolver="hogbom", gridder="standard",
                                   pbcor=True, savemodel="modelcolumn", imager_object=imager_objects[i], interactive=True))

        # Check minblperant, refant and spwmap
        selfcal_objects.append(Selfcal(visfile=clean_objects[i].inputvis,
                                       imagename=clean_objects[i].output, minblperant=2, refant="VA05", spwmap=[0, 0, 0, 0, 0, 0, 0, 0], Imager=clean_objects[i], want_plot=want_plot))

        phscal = Phasecal(minsnr=2.0, solint=solint_phs,
                          combine="spw", selfcal_object=selfcal_objects[i])

        phs_caltable = phscal.run()

    # ampcal=Ampcal(minsnr=2.0, solint=solint_amp, combine="scan",
    #                selfcal_object=parent_selfcal, input_caltable=phs_caltable)

    # amp_caltable=ampcal.run()

    # apcal = AmpPhasecal(minsnr=2.0, solint=solint_ap, combine="",
    #                    selfcal_object=parent_selfcal, input_caltable=amp_caltable)

    # apcal.run()

    # parent_selfcal.selfcal_output()
