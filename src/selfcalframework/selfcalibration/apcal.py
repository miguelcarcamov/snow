from casatasks import applycal, gaincal, rmtables

from .selfcal import Selfcal


class AmpPhasecal(Selfcal):

    def __init__(self, incremental: bool = False, solnorm: bool = True, **kwargs):
        super().__init__(**kwargs)

        self._calmode = 'ap'
        self._loops = len(self.solint)
        self.__incremental = incremental
        self.__solnorm = solnorm

        self._init_selfcal()

    def run(self):
        self._init_run("before_apcal")

        for i in range(0, self._loops):
            caltable = self.output_caltables + 'apcal_' + str(i)
            self._caltables.append(caltable)
            rmtables(caltable)

            self._set_attributes_from_dicts(i)

            if self.__incremental:
                gaincal(
                    vis=self.visfile,
                    field=self.imager.field,
                    caltable=caltable,
                    spw=self.imager.spw,
                    uvrange=self.uvrange,
                    gaintype=self.gaintype,
                    refant=self.refant,
                    calmode=self._calmode,
                    combine=self.combine,
                    solint=self.solint[i],
                    minsnr=self.minsnr,
                    minblperant=self.minblperant,
                    gaintable=self.input_caltable,
                    spwmap=self.spwmap,
                    solnorm=self.__solnorm
                )
            else:
                gaincal(
                    vis=self.visfile,
                    field=self.imager.field,
                    caltable=caltable,
                    spw=self.imager.spw,
                    uvrange=self.uvrange,
                    gaintype=self.gaintype,
                    refant=self.refant,
                    calmode=self._calmode,
                    combine=self.combine,
                    solint=self.solint[i],
                    minsnr=self.minsnr,
                    minblperant=self.minblperant,
                    spwmap=self.spwmap,
                    solnorm=self.__solnorm
                )

            self._plot_selfcal(
                caltable,
                xaxis="time",
                yaxis="amp",
                iteration="antenna",
                subplot=[4, 2],
                plotrange=[0, 0, 0.2, 1.8],
                want_plot=self.want_plot
            )

            version_name = 'before_apcal_' + str(i)
            self._save_selfcal(caltable_version=version_name, overwrite=True)
            self._caltables_versions.append(version_name)

            if self.__incremental:
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
            else:
                applycal(
                    vis=self.visfile,
                    spw=self.imager.spw,
                    spwmap=self.spwmap,
                    field=self.imager.field,
                    gaintable=[caltable],
                    gainfield='',
                    calwt=False,
                    flagbackup=False,
                    interp=self.interp,
                    applymode=self.applymode
                )
            if self.flag_dataset:
                self._flag_dataset(mode=self.flag_mode)

            self._run_imager(i)

            if self._finish_selfcal_iteration(i): break
