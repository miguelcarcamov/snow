from casatasks import tclean as tclean

from selfcalframework.utils.image_utils import *

from .imager import Imager


class Tclean(Imager):

    def __init__(
        self,
        nterms: int = 1,
        threshold: int = 0.0,
        nsigma: int = 0.0,
        interactive: bool = False,
        mask: str = "",
        use_mask: str = "auto-multithresh",
        negative_threshold: float = 0.0,
        low_noise_threshold: float = 1.5,
        noise_threshold: float = 4.25,
        sidelobe_threshold: float = 2.0,
        min_beam_frac: float = 0.3,
        specmode: str = "",
        gridder: str = "standard",
        wproj_planes: int = -1,
        deconvolver: str = "hogbom",
        uvtaper: list = [],
        scales: list = [],
        uvrange: str = "",
        pbcor: bool = False,
        cycle_niter: int = 0,
        clean_savemodel: str = None,
        **kwargs
    ):
        super().__init__(**kwargs)
        self.name = "TClean"
        self.nterms = nterms
        self.threshold = threshold
        self.nsigma = nsigma
        self.interactive = interactive
        self.mask = mask
        self.use_mask = use_mask
        self.negative_threshold = negative_threshold
        self.low_noise_threshold = low_noise_threshold
        self.noise_threshold = noise_threshold
        self.sidelobe_threshold = sidelobe_threshold
        self.min_beam_frac = min_beam_frac
        self.specmode = specmode
        self.gridder = gridder
        self.wproj_planes = wproj_planes
        self.deconvolver = deconvolver
        self.uvtaper = uvtaper
        self.scales = scales
        self.uvrange = uvrange
        self.pbcor = pbcor
        self.cycle_niter = cycle_niter
        self.clean_savemodel = clean_savemodel

        if self.save_model:
            self.clean_savemodel = "modelcolumn"

    def run(self, imagename=""):
        __imsize = [self.M, self.N]
        tclean(
            vis=self.inputvis,
            imagename=imagename,
            field=self.field,
            phasecenter=self.phase_center,
            uvrange=self.uvrange,
            datacolumn=self.data_column,
            specmode=self.specmode,
            stokes=self.stokes,
            deconvolver=self.deconvolver,
            scales=self.scales,
            nterms=self.nterms,
            imsize=__imsize,
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
            negativethreshold=self.negative_threshold,
            lownoisethreshold=self.low_noise_threshold,
            noisethreshold=self.noise_threshold,
            sidelobethreshold=self.sidelobe_threshold,
            minbeamfrac=self.min_beam_frac,
            cycleniter=self.cycle_niter,
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
