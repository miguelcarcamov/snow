import os
import shlex
import subprocess

from casatasks import exportfits, fixvis, imhead, immath, importfits, tclean
from casatools import image, quanta

from .imager import Imager


class GPUvmem(Imager):

    def __init__(
        self,
        executable="gpuvmem",
        gpublocks=[16, 16, 256],
        initialvalues=[],
        regfactors=[],
        gpuids=[0],
        residualoutput="residuals.ms",
        model_input="",
        modelout="mod_out.fits",
        user_mask="",
        force_noise=None,
        griddingthreads=4,
        positivity=True,
        ftol=1e-12,
        noise_cut=10.0,
        gridding=False,
        printimages=False,
        **kwargs
    ):
        super(GPUvmem, self).__init__(**kwargs)
        self.name = "GPUvmem"
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        # self.__dict__.update(kwargs)
        if self.phasecenter != "":
            fixvis(
                vis=self.visfile,
                outputvis=self.visfile,
                field=self.field,
                phasecenter=self.phasecenter
            )

    def getRegfactors(self):
        return self.regfactors

    def setRegFactors(self, regfactors):
        self.factors = regfactors

    def _restore(self, model_fits="", residual_ms="", restored_image="restored"):
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

    def _make_canvas(self, name="model_input"):
        fitsimage = name + '.fits'
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
        exportfits(imagename=name + '.image', fitsimage=fitsimage, overwrite=True)
        return fitsimage

    def run(self, imagename=""):
        if self.model_input == "":
            self.model_input = self._make_canvas(imagename + "_input")
        model_output = imagename + ".fits"
        residual_output = imagename + "_" + self.residualoutput
        restored_image = imagename + ".restored"

        args = self.executable + " -X " + str(self.gpublocks[0]) + " -Y " + str(self.gpublocks[1]) + " -V " + str(
            self.gpublocks[2]) \
               + " -i " + self.inputvis + " -o " + residual_output + " -z " + ",".join(map(str, self.initialvalues)) \
               + " -Z " + ",".join(map(str, self.regfactors)) + " -G " + ",".join(map(str, self.gpuids)) \
               + " -m " + self.model_input + " -O " + model_output + " -N " + str(self.noise_cut) \
               + " -R " + str(self.robust) + " -t " + str(self.niter)

        if self.user_mask != "":
            args += " -U " + self.user_mask

        if self.force_noise is not None:
            args += " -n " + str(self.force_noise)

        if self.gridding:
            args += " -g " + str(self.griddingthreads)

        if self.printimages:
            args += " --print-images"

        if not self.positivity:
            args += " --nopositivity"

        if self.verbose:
            args += " --verbose"

        if self.savemodel:
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
            residual_fits, restored_fits = self._restore(
                model_fits=model_output, residual_ms=residual_output, restored_image=restored_image
            )

        # Calculate SNR and standard deviation
        self.calculateStatistics_FITS(
            signal_fits_name=restored_fits, residual_fits_name=residual_fits
        )
