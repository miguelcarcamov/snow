import numpy as np
import imager
import selfcal
import sys

visfile = sys.argv[3]
output = sys.argv[4]

clean_imager = Clean(specmode="mfs", deconvolver="hogbom", gridder="standard", pbcor=True, savemodel="modelcolumn", inputvis=visfile, output=output,
                     niter=100, M=1024, N=1024, deltax="0.02arcsec", stokes="I", datacolumn="corrected", robust=0.5)
parent_selfcal = Selfcal(visfile=clean_imager.inputvis,
                         imagename=clean_imager.output, minblperant=3, refant="Kn,Mk2,Cm,Pi,De,Da", spwmap=[0, 0, 0, 0, 0, 0, 0, 0], clean_imager)

solint_phs = ['inf', '30.25s', 'int']
solint_ap = ['inf']

phscal = Phasecal(minsnr=2.0, solint=solint_phs, combine="spw", parent_selfcal)
ampcal = Ampcal(minsnr=2.0, solint=solint_ap, combine="spw")
