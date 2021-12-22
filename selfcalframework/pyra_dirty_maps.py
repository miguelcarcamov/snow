import sys
from astropy.io import fits
from pyralysis.io import DaskMS
# from pyralysis.reconstruction import Image
import astropy.units as u
import numpy as np
from pyralysis.units import lambdas_equivalencies
from pyralysis.transformers.weighting_schemes import Uniform, Robust
from pyralysis.transformers import Gridder, HermitianSymmetry, DirtyMapper
from pyralysis.io import FITS
from astropy.units import Quantity
import re

#python pyra_dirty_maps.py residual_ms file_pyra_dirty weighting robust side cell

residual_ms=sys.argv[1]
file_pyra_dirty=sys.argv[2]
weighting=sys.argv[3] # string
robust=float(sys.argv[4])
side=int(sys.argv[5])
pixscl=float(sys.argv[6]) # arcsec
model_image=sys.argv[7]

print("pyra: loading residual ms")
x = DaskMS(input_name=residual_ms)
dataset = x.read()
#dataset.field.mean_ref_dir
#dataset.psf[0].sigma

print("pyra: done")

# h_symmetry = HermitianSymmetry(input_data=dataset)
# h_symmetry.apply()

imsize = [side, side]

#dx = Quantity([dataset.theo_resolution/7, dataset.theo_resolution/7])
#dx.to(u.arcsec)
dx = [pixscl, pixscl] * u.arcsec

du = (1/(imsize*dx)).to(u.lambdas, equivalencies=lambdas_equivalencies())

gridder = Gridder(imsize=imsize, cellsize=dx, hermitian_symmetry=True)

fits_io = FITS()

dirty_mapper = DirtyMapper(input_data=dataset, imsize=imsize, cellsize=dx, stokes="I", hermitian_symmetry=False)

dirty_images_natural = dirty_mapper.transform()
dirty_image_natural = dirty_images_natural[0].data[0].compute()
#dirty_beam_natural = dirty_images_natural[1].data[0].compute()

fits_io.write(dirty_images_natural[0].data, output_name='pyra_dirty_nat.fits')

print("selecting weighting schemes")
if re.search('Natural',weighting,re.IGNORECASE):
    fits_io.write(dirty_images_natural[0].data, output_name=file_pyra_dirty)

    #gridded_visibilities_nat = dirty_mapper.uvgridded_visibilities.compute()
    #gridded_weights_nat = dirty_mapper.uvgridded_weights.compute()

elif re.search('Briggs',weighting,re.IGNORECASE):
    ######################################################################
    # Briggs weighting:

    robust = Robust(input_data=dataset, robust_parameter=robust, gridder=gridder)
    robust.apply()

    dataset.calculate_psf()
    
    dirty_mapper = DirtyMapper(input_data=dataset, imsize=imsize, cellsize=dx, stokes="I", hermitian_symmetry=False)

    dirty_images_robust = dirty_mapper.transform()
    dirty_image = dirty_images_robust[0].data[0].compute()
    #dirty_beam = dirty_images_robust[1].data[0].compute()
    fits_io.write(dirty_images_robust[0].data, output_name=file_pyra_dirty)


elif re.search('Uniform',weighting,re.IGNORECASE):
    ######################################################################
    # Uniform weighting

    uniform = Uniform(input_data=dataset, gridder=Gridder(imsize=imsize, uvcellsize=du))

    uniform.apply()

    dataset.calculate_psf()
    
    dirty_mapper = DirtyMapper(input_data=dataset, imsize=imsize, cellsize=dx, stokes="I", hermitian_symmetry=False)
    
    dirty_images_uniform = dirty_mapper.transform()
    
    dirty_image = dirty_images_uniform[0].data[0].compute()
    dirty_beam = dirty_images_uniform[1].data[0].compute()
    
    fits_io.write(dirty_images_uniform[0].data, output_name=file_pyra_dirty)

else:
    print("unrecognised weighting scheme:",weighting)


print("pyra: pasting input header")
print("input image from file ",file_pyra_dirty)
hdu_residuals=fits.open(file_pyra_dirty)
hdr_res=hdu_residuals[0].header
hdu_model=fits.open(model_image)
print("using header from file ",model_image)
hdr=hdu_model[0].header
hdr['BUNIT']='Jy/beam'
hdr['BMAJ']=hdr_res['BMAJ']
hdr['BMIN']=hdr_res['BMIN']
hdr['BPA']=hdr_res['BPA']

hdu_residuals[0].header=hdr

hdu_residuals.writeto(file_pyra_dirty,overwrite=True)
print("pyra: exit")
