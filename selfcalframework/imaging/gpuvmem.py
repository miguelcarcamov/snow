import os
import shlex
import subprocess
from pathlib import Path

import numpy as np
from casatasks import exportfits, fixvis, imhead, immath, importfits, tclean
from casatools import image, quanta
from reproject import reproject_interp

from ..utils.image_utils import get_hdul, get_header
from .imager import Imager


class GPUvmem(Imager):

    def __init__(
        self,
        executable: str = "gpuvmem",
        gpu_blocks: list = [16, 16, 256],
        initial_values: list = [],
        regfactors: list = [],
        gpuids: list = [0],
        residual_output: str = "residuals.ms",
        model_input: str = None,
        model_out: str = "mod_out.fits",
        user_mask: str = None,
        force_noise: float = None,
        gridding_threads: int = 4,
        positivity: bool = True,
        ftol: float = 1e-12,
        noise_cut: float = 10.0,
        gridding: bool = False,
        print_images: bool = False,
        **kwargs
    ):
        super().__init__(**kwargs)
        self.name = "GPUvmem"
        self.executable = executable
        self.gpu_blocks = gpu_blocks
        self.initial_values = initial_values
        self.regfactors = regfactors
        self.gpuids = gpuids
        self.residual_output = residual_output
        self.model_input = model_input
        self.model_out = model_out
        self.user_mask = user_mask
        self.force_noise = force_noise
        self.gridding_threads = gridding_threads
        self.positivity = positivity
        self.ftol = ftol
        self.noise_cut = noise_cut
        self.gridding = gridding
        self.print_images = print_images

        if self.phase_center != "":
            fixvis(
                vis=self.visfile,
                outputvis=self.visfile,
                field=self.field,
                phasecenter=self.phase_center
            )

    @property
    def user_mask(self):
        return self.__user_mask

    @user_mask.setter
    def user_mask(self, user_mask):
        self.__user_mask = user_mask
        if self.__user_mask == "":
            pass
        else:
            self.__check_mask()

    def __check_mask(self, order="bilinear"):

        if self.user_mask is not None and self.model_input is not None:
            if os.path.exists(self.user_mask) and os.path.exists(self.model_input):
                header_mask = get_header(self.user_mask)
                header_model = get_header(self.model_input)

                if np.fabs(
                    (header_mask['CDELT2'] - header_model['CDELT2']) / header_model['CDELT2']
                ) < 1E-3 or header_mask['NAXIS1'] != header_model['NAXIS1']:
                    print("The mask header is not the same as the model image, resampling...")
                    array, footprint = reproject_interp(self.user_mask, header_model, order=order)
                    path_object = Path(self.user_mask)
                    resampled_mask_name = "{0}_{2}{1}".format(
                        Path.joinpath(path_object.parent, path_object.stem), path_object.suffix,
                        "resampled"
                    )
                    fits.writeto(resampled_mask_name, array, header_model, overwrite=True)
                    self.user_mask = resampled_mask_name
            else:
                print("User mask: {0}".format(self.user_mask))
                print("Model input: {0}".format(self.model_input))
                raise FileNotFoundError("Either user mask of model input does not exist")

    def __restore(self, model_fits="", residual_ms="", restored_image="restored"):
        qa = quanta()
        ia = image()
        residual_image = residual_ms.partition(".ms")[0] + ".residual"
        os.system(
            "rm -rf *.log *.last " + residual_image +
            ".* mod_out convolved_mod_out convolved_mod_out.fits " + restored_image + " " +
            restored_image + ".fits"
        )

        importfits(imagename="model_out", fitsimage=model_fits, overwrite=True)
        shape = imhead(imagename="model_out", mode="get", hdkey="shape")
        pix_num = shape[0]
        cdelt = imhead(imagename="model_out", mode="get", hdkey="cdelt2")
        cdelta = qa.convert(v=cdelt, outunit="arcsec")
        cdeltd = qa.convert(v=cdelt, outunit="deg")
        pix_size = str(cdelta['value']) + "arcsec"

        tclean(
            vis=residual_ms,
            imagename=residual_image,
            specmode='mfs',
            deconvolver='hogbom',
            niter=0,
            stokes=self.stokes,
            nterms=1,
            weighting=self.weighting,
            robust=self.robust,
            imsize=[self.M, self.N],
            cell=self.cell,
            datacolumn='data'
        )

        exportfits(
            imagename=residual_image + ".image",
            fitsimage=residual_image + ".image.fits",
            overwrite=True,
            history=False
        )

        ia.open(infile=residual_image + ".image")
        rbeam = ia.restoringbeam()
        ia.done()
        ia.close()

        bmaj = imhead(imagename=residual_image + ".image", mode="get", hdkey="beammajor")
        bmin = imhead(imagename=residual_image + ".image", mode="get", hdkey="beamminor")
        bpa = imhead(imagename=residual_image + ".image", mode="get", hdkey="beampa")

        minor = qa.convert(v=bmin, outunit="deg")
        pa = qa.convert(v=bpa, outunit="deg")

        ia.open(infile="model_out")
        im2 = ia.convolve2d(
            outfile="convolved_model_out",
            axes=[0, 1],
            type='gauss',
            major=bmaj,
            minor=bmin,
            pa=bpa,
            overwrite=True
        )
        im2.done()
        ia.done()
        ia.close()

        exportfits(
            imagename="convolved_model_out",
            fitsimage="convolved_model_out.fits",
            overwrite=True,
            history=False
        )
        ia.open(infile="convolved_model_out.fits")
        ia.setrestoringbeam(beam=rbeam)
        ia.done()
        ia.close()

        imagearr = ["convolved_model_out.fits", residual_image + ".image.fits"]

        immath(imagename=imagearr, expr=" (IM0   + IM1) ", outfile=restored_image)

        exportfits(
            imagename=restored_image,
            fitsimage=restored_image + ".fits",
            overwrite=True,
            history=False
        )

        return residual_image + ".image.fits", restored_image + ".fits"

    def __make_canvas(self, name="model_input"):
        fits_image = name + '.fits'
        tclean(
            vis=self.inputvis,
            imagename=name,
            specmode='mfs',
            niter=0,
            deconvolver='hogbom',
            interactive=False,
            cell=self.cell,
            stokes=self.stokes,
            robust=self.robust,
            imsize=[self.M, self.N],
            weighting=self.weighting
        )
        exportfits(imagename=name + '.image', fitsimage=fits_image, overwrite=True)
        return fits_image

    def run(self, imagename=""):
        if self.model_input is None:
            self.model_input = self.__make_canvas(imagename + "_input")
        model_output = imagename + ".fits"
        _residual_output = imagename + "_" + self.residual_output
        restored_image = imagename + ".restored"

        args = self.executable + " -X " + str(self.gpu_blocks[0]) + " -Y " + str(self.gpu_blocks[1]) + " -V " + str(
            self.gpu_blocks[2]) \
               + " -i " + self.inputvis + " -o " + _residual_output + " -z " + ",".join(map(str, self.initial_values)) \
               + " -Z " + ",".join(map(str, self.regfactors)) + " -G " + ",".join(map(str, self.gpuids)) \
               + " -m " + self.model_input + " -O " + model_output + " -N " + str(self.noise_cut) \
               + " -R " + str(self.robust) + " -t " + str(self.niter)

        if self.user_mask is not None and type(self.user_mask) is str:
            args += " -U " + self.user_mask

        if self.force_noise is not None:
            args += " -n " + str(self.force_noise)

        if self.gridding:
            args += " -g " + str(self.gridding_threads)

        if self.print_images:
            args += " --print-images"

        if not self.positivity:
            args += " --nopositivity"

        if self.verbose:
            args += " --verbose"

        if self.save_model:
            args += " --save_modelcolumn"

        print(args)
        args = shlex.split(args)
        print(args)

        # Run gpuvmem and wait until it finishes
        p = subprocess.Popen(args, env=os.environ)
        p.wait()

        if not os.path.exists(model_output):
            raise FileNotFoundError("The model image has not been created")
        else:
            # Restore the image
            residual_fits, restored_fits = self.__restore(
                model_fits=model_output,
                residual_ms=_residual_output,
                restored_image=restored_image
            )

        # Calculate SNR and standard deviation
        self._calculate_statistics_fits(
            signal_fits_name=restored_fits, residual_fits_name=residual_fits
        )
