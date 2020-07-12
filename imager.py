import numpy as np
from tclean import tclean


class Imager(object):

    def __init__(self, inputvis="", output="", niter=0, M=0, N=0, deltax="", stokes="", datacolumn="", robust=0, field="", spw=""):
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
        self.spw = spw

    def getField(self):
        return self.field

    def getSpw(self):
        return self.spw


class Clean(Imager):
    def __init__(self, specmode="", deconvolver="", nterms=1, threshold=0.0, interactive=False, gridder="", pbcor=False,
                 savemodel="", usepointing=False, usemask="", negativethreshold=0.0, lownoisethreshold=1.5 , noisethreshold=4.25 ,
                 sidelobethreshold=2.0, minbeamfrac=0.3 , imager_object=None):
        super(Clean, self).__init__(imager_object.inputvis, imager_object.output, imager_object.niter, imager_object.M, imager_object.N,
                                    imager_object.deltax, imager_object.stokes, imager_object.datacolumn, imager_object.robust, imager_object.field)
        self.specmode = specmode
        self.deconvolver = deconvolver
        self.nterms = nterms
        self.threshold = threshold
        self.interactive = interactive
        self.gridder = gridder
        self.pbcor = pbcor
        self.savemodel = savemodel
        self.usepointing = usepointing
        self.usemask = usemask
        self.negativethreshold = negativethreshold
        self.lownoisethreshold = lownoisethreshold
        self.noisethreshold = noisethreshold
        self.sidelobethreshold = sidelobethreshold
        self.minbeamfrac = minbeamfrac

    def run(self, imagename=""):
        imsize = [self.M, self.N]
        tclean(vis=self.inputvis, imagename=imagename, field=self.field,
               datacolumn=self.datacolumn, specmode=self.specmode, stokes=self.stokes, deconvolver=self.deconvolver, nterms=self.nterms,
               imsize=imsize, cell=self.deltax, weighting="briggs", robust=self.robust, niter=self.niter, threshold=self.threshold,
               interactive=self.interactive, gridder=self.gridder, pbcor=self.pbcor, savemodel=self.savemodel, usepointing=self.usepointing,
               usemask=self.usemask, negativethreshold=self.negativethreshold, lownoisethreshold=self.lownoisethreshold, noisethreshold=self.noisethreshold,
               sidelobethreshold=self.sidelobethreshold, minbeamfrac=self.minbeamfrac)


class WSClean(Imager):
    def __init__(self, imager_object=None):
        super(Clean, self).__init__(imager_object.inputvis, imager_object.output, imager_object.niter, imager_object.M, imager_object.N,
                                    imager_object.deltax, imager_object.stokes, imager_object.datacolumn, imager_object.robust, imager_object.field)
