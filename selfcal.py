import numpy as np

class Selfcal:

    def __init__(self, visfile="", imagename="", minblperant=0, refant-"", spwmap=[], Imager=None):
        self.visfile=visfile
        self.imagename=imagename
        self.minblperant = minblperant
        self.refant = refant
        self.spwmap = spwmap
        self.Imager = Imager

class Ampcal(Selfcal):
    def __init__(self, minsnr=1.0, solint=[], combine="", phs_caltable="", visfile="", imagename="", minblperant=0, refant="", spwmap=[], Imager=None):
        super().__init__(visfile, imagename, minblperant, refant, spwmap, Imager)
        self.calmode = 'ap'
        self.minsnr = minsnr
        self.solint = solint
        self.combine = combine
        self.phs_caltable = phs_caltable
        self.loops = len(self.solint)

    def run(self):
        for i in range(0, self.loops):
            caltable = 'apcal_'+str(i)
            rmtables(caltable)
            gaincal(vis=self.visfile, field=Imager.getField(), caltable = caltable, spw=Imager.getspw(), gaintype ='G', refant = self.refant, calmode = self.calmode, combine = self.combine, solint=self.solint[i], minsnr=self.minsnr, minblperant=minblperant, gaintable=self.phs_caltable, spwmap=self.spwmap, solnorm=True)

            plot_selfcal(caltable)

            applycal(vis=self.visfile, spwmap=[self.spwmap, self.spwmap], field=Imager.getField(), gaintable=[self.phs_caltable, caltable], gainfield='', calwt=False, flagbackup=False, interp='linearperobs')

            flagmanager(vis=self.visfile, mode='save', versionname='after_apcal'+str(i))

            imagename = imagename + '_ap'+str(i)

            Imager.run(imagename)


class Phasecal(Selfcal):
    def __init__(self, minsnr=1.0, solint=[], combine="", visfile="", imagename="", minblperant=0, refant="", spwmap=[], Imager=None):
        super().__init__(visfile, imagename, minblperant, refant, spwmap, Imager)
        self.calmode = 'p'
        self.minsnr = minsnr
        self.solint = solint
        self.combine = combine
        self.loops = len(self.solint)

    def run(self):
        flagmanager(vis=self.visfile, mode='save', versionname='before_phasecal', merge='replace')
        caltable = ""
        for i in range(0, self.loops):
            imagename = self.imagename + '_ph' + str(i)
            Imager.run(imagename)
            caltable = 'pcal'+str(i)
            rmtables(caltable)

            gaincal(vis=self.visfile, caltable=caltable, field=Imager.getField(), spw=Imager.getspw(), gaintype='G', refant=self.refant, calmode=self.calmode, combine=self.combine, solint=self.solint[i], minsnr=self.minsnr, minblperant=self.minblperant)

            plot_selfcal(caltable)

            applycal(vis=self.visfile, field=Imager.getField(), spwmap=self.spwmap, gaintable=[caltable], gainfield='', calwt=False, flagbackup=False, interp='linearperobs')

            flagmanager(vis=self.visfile, mode='save', versionname='after_'+caltable)
        return caltable


