import numpy as np


class Imager:

    def __init__(self, inputvis="", output="", niter=0, M=0, N=0, deltax="", stokes="", datacolumn="", robust=0, field=""):
        self.inputvis = inputvis
        self.output = output
        self.niter = niter
        self.M = M
        self.N = N
        self.deltax = deltax
        self.stokes = stokes
        self.datacolumn = datacolumn
        self.robust = robust
        self.field = field

    def getField(self):
        return self.field


class Clean(Imager):
    def __init__(self, specmode="", deconvolver="", nterms=1, threshold=0.0, interactive=False, gridder="", pbcor=False,
                 savemodel="", usepointing=False, inputvis="", output="", niter=0, M=0, N=0, deltax="", stokes="", datacolumn="", robust=0, field=""):
        super().__init__(inputvis, output, niter, M, N,
                         deltax, stokes, datacolumn, robust, field)
        self.specmode = specmode
        self.deconvolver = deconvolver
        self.nterms = nterms
        self.threshold = threshold
        self.interactive = interactive
        self.gridder = gridder
        self.pbcor = pbcor
        self.savemodel = savemodel
        self.usepointing = usepointing

    def run(self, imagename=""):
        imsize = [self.M, self.N]
        tclean(vis=self.inputvis, imagename=imagename, field=self.field,
               datacolumn=self.datacolumn, specmode=self.specmode, stokes=self.stokes, deconvolver=self.deconvolver, nterms=self.nterms,
               imsize=imsize, cell=self.deltax, weighting="briggs", robust=self.robust, niter=self.niter, threshold=self.threshold,
               interactive=self.interactive, gridder=self.gridder, pbcor=self.pbcor, savemodel=self.savemodel, usepointing=self.usepointing)
