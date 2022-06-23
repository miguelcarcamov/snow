import os
import shutil
from pathlib import Path

from casatasks import applycal, gaincal, rmtables

from .selfcal import Selfcal


class Ampcal(Selfcal):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self._calmode = 'a'
        self._loops = len(self.solint)

        self._init_selfcal()

    def run(self):
        caltable = ""
        self._init_run("before_ampcal")

        for i in range(0, self._loops):
            caltable = self.output_caltables + 'ampcal_' + str(i)
            self._caltables.append(caltable)
            rmtables(caltable)

            self._set_attributes_from_dicts(i)

            gaincal(
                vis=self.visfile,
                field=self.imager.field,
                caltable=caltable,
                spw=self.imager.spw,
                uvrange=self.uvrange,
                gaintype=self.gaintype,
                refant=self.refant,
                calmode=self.calmode,
                combine=self.combine,
                solint=self.solint[i],
                minsnr=self.minsnr,
                minblperant=self.minblperant,
                gaintable=self.input_caltable,
                spwmap=self.spwmap,
                solnorm=True
            )

            self.plot_selfcal(
                caltable,
                xaxis="time",
                yaxis="amp",
                iteration="antenna",
                subplot=[4, 2],
                plotrange=[0, 0, 0.2, 1.8],
                want_plot=self.want_plot
            )

            versionname = 'before_ampcal_' + str(i)
            self._save_selfcal(caltable_version=versionname, overwrite=True)
            self._caltables_versions.append(versionname)

            applycal(
                vis=self.visfile,
                spw=self.imager.spw,
                spwmap=[self.spwmap, self.spwmap],
                field=self.imager.field,
                gaintable=[self.input_caltable, caltable],
                gainfield='',
                calwt=False,
                flagbackup=False,
                interp=self.interp,
                applymode=self.applymode
            )
            self.input_caltable = caltable

            if self.flag_dataset_bool:
                self._flag_dataset(mode=self.flag_mode)

            self._run_imager(i)

            if self._finish_selfcal_loop(i): break
