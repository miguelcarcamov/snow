import numpy as np
from tclean import tclean


class Imager(object):

    def __init__(self, robust=2.0, field="", spw="", **kwargs):
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(kwargs)

    def getField(self):
        return self.field

    def getSpw(self):
        return self.spw

    def getVis(self):
        return self.inputvis

    def getOutputPath(self):
        return self.output


class Clean(Imager):
    def __init__(self, nterms=1, threshold=0.0, interactive=False, usemask="auto-multithresh", negativethreshold=0.0, lownoisethreshold=1.5, noisethreshold=4.25,
                 sidelobethreshold=2.0, minbeamfrac=0.3, **kwargs):
        super(Clean, self).__init__(**kwargs)
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(kwargs)

    def run(self, imagename=""):
        imsize = [self.M, self.N]
        tclean(vis=self.inputvis, imagename=imagename, field=self.field,
               datacolumn=self.datacolumn, specmode=self.specmode, stokes=self.stokes, deconvolver=self.deconvolver, nterms=self.nterms,
               imsize=imsize, cell=self.cell, weighting="briggs", robust=self.robust, niter=self.niter, threshold=self.threshold,
               interactive=self.interactive, gridder=self.gridder, pbcor=self.pbcor, savemodel=self.savemodel, usemask=self.usemask,
               negativethreshold=self.negativethreshold, lownoisethreshold=self.lownoisethreshold, noisethreshold=self.noisethreshold,
               sidelobethreshold=self.sidelobethreshold, minbeamfrac=self.minbeamfrac)


class WSClean(Imager):
    def __init__(self, **kwargs):
        super(Clean, self).__init__(**kwargs)
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(kwargs)
