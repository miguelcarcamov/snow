import os
import shutil

from casatasks import applycal, gaincal, rmtables

from .selfcal import Selfcal


class Phasecal(Selfcal):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self._calmode = 'p'
        self._loops = len(self.solint)

        self._init_selfcal()

    def run(self):
        caltable = "before_selfcal"
        self._save_selfcal(caltable_version=caltable, overwrite=True)
        self.caltables_versions.append(caltable)
        self._init_run("_original")

        for i in range(0, self._loops):
            caltable = 'pcal' + str(i)
            self._caltables.append(caltable)
            rmtables(caltable)

            self._set_attributes_from_dicts(i)

            gaincal(
                vis=self.visfile,
                caltable=caltable,
                field=self.imager.field,
                spw=self.imager.spw,
                uvrange=self.uvrange,
                gaintype=self.gaintype,
                refant=self.refant,
                calmode=self._calmode,
                combine=self.combine,
                solint=self.solint[i],
                minsnr=self.minsnr,
                spwmap=self.spwmap,
                minblperant=self.minblperant
            )

            self._plot_selfcal(
                caltable,
                xaxis="time",
                yaxis="phase",
                iteration="antenna",
                subplot=[4, 2],
                plotrange=[0, 0, -180, 180],
                want_plot=self.want_plot
            )

            versionname = 'before_phasecal_' + str(i)
            self._save_selfcal(caltable_version=versionname, overwrite=True)
            self._caltables_versions.append(versionname)

            applycal(
                vis=self.visfile,
                field=self.field,
                spw=self.spw,
                spwmap=self.spwmap,
                gaintable=[caltable],
                gainfield='',
                calwt=False,
                flagbackup=False,
                interp=self.interp,
                applymode=self.applymode
            )

            if self.flag_dataset:
                self._flag_dataset(mode=self.flag_mode)

            imagename = self.image_name + '_ph' + str(i)

            self.imager.run(imagename)

            self._psnr_history.append(self.imager.psnr)

            print("Solint: " + str(self.solint[i]) + " - PSNR: " + str(self._psnr_history[-1]))
            print("Noise: " + str(self.imager.stdv * 1000.0) + " mJy/beam")

            if not self._finish_selfcal_loop(i): break
