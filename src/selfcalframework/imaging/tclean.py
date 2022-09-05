from casatasks import tclean

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
        """
        tclean imager object

        Parameters
        ----------
        nterms :
            Number of Taylor coefficients in the spectral model
        threshold :
            Stopping threshold (number in units of Jy, or string)
        nsigma :
            Multiplicative factor for rms-based threshold stopping. N-sigma threshold is calculated as
            nsigma * rms value per image plane determined from a robust statistics
        interactive :
            Modify masks and parameters at runtime
        mask :
            Mask (a list of image name(s) or region file(s) or region string(s) .The name of a CASA image
            or region file or region string that specifies a 1/0 mask to be used for deconvolution.
        use_mask :
            Type of mask(s) to be used for deconvolution
        negative_threshold :
            Sub-parameter for “auto-multithresh”: mask threshold for negative
            features: -1.0* negativethreshold * rms + location(=median)
        low_noise_threshold :
            Sub-parameter for “auto-multithresh”: mask threshold to grow previously masked regions
            via binary dilation: lownoisethreshold * rms in residual image + location (=median)
        noise_threshold :
            Sub-parameter for “auto-multithresh”: mask threshold based on the noise
            level: noisethreshold * rms + location (=median)
        sidelobe_threshold :
            Sub-parameter for “auto-multithresh”: mask threshold based on sidelobe
            levels: sidelobethreshold * max_sidelobe_level * peak residual
        min_beam_frac :
            Sub-parameter for “auto-multithresh”: minimum beam fraction in size to prune masks smaller than
            minbeamfrac * beam<=0.0 : No pruning
        specmode :
            Spectral definition mode (mfs, cube, cubedata)
        gridder :
            Gridding options (standard, wproject, widefield, mosaic, awproject)
        wproj_planes :
            Number of distinct w-values at which to compute and use different gridding convolution
            functions for W-Projection
        deconvolver :
            Name of minor cycle algorithm (hogbom, clark, multiscale, mtmfs, mem, clarkstokes)
        uvtaper :
            uv-taper on outer baselines in uv-plane
        scales :
            List of scale sizes (in pixels) for multi-scale and mtmfs algorithms.
        uvrange :
            Select data within uvrange (default unit is meters)
        pbcor :
            Apply PB correction on the output restored image
        cycle_niter :
            Maximum number of minor-cycle iterations (per plane) before triggering a major cycle
        clean_savemodel :
            Options to save model visibilities (none, virtual, modelcolumn)
        kwargs :
            General imager arguments
        """
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

        self._calculate_statistics_msimage(
            signal_ms_name=restored_image, residual_ms_name=residual_image
        )
