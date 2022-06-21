from casatasks import tclean

from selfcalframework.utils.image_utils import *

from .imager import Imager


class Clean(Imager):

    def __init__(
        self,
        nterms=1,
        threshold=0.0,
        nsigma=0.0,
        interactive=False,
        mask="",
        usemask="auto-multithresh",
        negativethreshold=0.0,
        lownoisethreshold=1.5,
        noisethreshold=4.25,
        sidelobethreshold=2.0,
        minbeamfrac=0.3,
        specmode="",
        gridder="standard",
        wprojplanes=-1,
        deconvolver="hogbom",
        uvtaper=[],
        scales=[],
        uvrange="",
        pbcor=False,
        cycleniter=0,
        clean_savemodel=None,
        **kwargs
    ):
        super(Clean, self).__init__(**kwargs)
        self.name = "TClean"
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        # self.__dict__.update(kwargs)

        if (self.savemodel):
            self.clean_savemodel = "modelcolumn"

    def getNSigma(self):
        return self.nsigma

    def setNSigma(self, nsigma):
        self.nsigma = nsigma

    def getThreshold(self):
        return self.threshold

    def setThreshold(self):
        return self.threshold

    def run(self, imagename=""):
        imsize = [self.M, self.N]
        tclean(
            vis=self.inputvis,
            imagename=imagename,
            field=self.field,
            phasecenter=self.phasecenter,
            uvrange=self.uvrange,
            datacolumn=self.datacolumn,
            specmode=self.specmode,
            stokes=self.stokes,
            deconvolver=self.deconvolver,
            scales=self.scales,
            nterms=self.nterms,
            imsize=imsize,
            cell=self.cell,
            weighting=self.weighting,
            robust=self.robust,
            niter=self.niter,
            threshold=self.threshold,
            nsigma=self.nsigma,
            interactive=self.interactive,
            gridder=self.gridder,
            mask=self.mask,
            pbcor=self.pbcor,
            uvtaper=self.uvtaper,
            savemodel=self.clean_savemodel,
            usemask=self.usemask,
            negativethreshold=self.negativethreshold,
            lownoisethreshold=self.lownoisethreshold,
            noisethreshold=self.noisethreshold,
            sidelobethreshold=self.sidelobethreshold,
            minbeamfrac=self.minbeamfrac,
            cycleniter=self.cycleniter,
            verbose=self.verbose
        )

        if self.deconvolver != "mtmfs":
            restored_image = imagename + ".image"
            residual_image = imagename + ".residual"
        else:
            restored_image = imagename + ".image.tt0"
            residual_image = imagename + ".residual.tt0"

        self.calculateStatistics_MSImage(
            signal_ms_name=restored_image, residual_ms_name=residual_image
        )
