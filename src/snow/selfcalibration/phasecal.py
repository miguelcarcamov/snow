from dataclasses import dataclass

from casatasks import applycal, gaincal, rmtables

from .selfcal import Selfcal


@dataclass(init=False, repr=True)
class Phasecal(Selfcal):

    def __init__(self, **kwargs):
        """
        Phase only self-calibration object

        Parameters
        ----------
        kwargs :
            General self-calibration arguments
        """
        super().__init__(**kwargs)

        self._calmode = 'p'
        self._loops = len(self.solint)

        self._copy_directory_at_start()
        self._init_selfcal()

    def run(self):
        caltable = "before_selfcal"
        self._save_selfcal(caltable_version=caltable, overwrite=True)
        self._caltables_versions.append(caltable)
        self._init_run("_original")

        for i in range(0, self._loops):
            caltable = self.output_caltables + 'pcal' + str(i)
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

            version_name = 'before_phasecal_' + str(i)
            self._save_selfcal(caltable_version=version_name, overwrite=True)
            self._caltables_versions.append(version_name)

            print("Applying calibration tables to {0} file".format(self.visfile))
            applycal(
                vis=self.visfile,
                field=self.imager.field,
                spw=self.imager.spw,
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

            self._run_imager(i)

            if self._finish_selfcal_iteration(i):
                break
