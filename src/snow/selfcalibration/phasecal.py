from dataclasses import dataclass

from casatasks import rmtables

from .caltable import CalStack, CalTable
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

        # Phasecal solves and applies only its own tables, so earlier amplitude solutions or an
        # input table would be silently discarded (#25).
        if self.input_caltable:
            raise ValueError(
                "Phasecal does not pre-apply input_caltable ({0}); its solutions would be "
                "discarded.".format(self.input_caltable)
            )
        previous = self.previous_selfcal
        if previous is not None and any(
            caltable.solves_amplitude for caltable in previous._caltable_stack
        ):
            raise ValueError(
                "Phasecal after an amplitude self-calibration ({0}) would discard its amplitude "
                "solutions.".format(type(previous).__name__)
            )

        self._copy_directory_at_start()
        self._init_selfcal()

    def _preapplied_caltables(self) -> CalStack:
        """
        Phase self-calibration solves from scratch against the model: earlier phase tables are
        not pre-applied but replaced
        """
        return CalStack()

    def _accept_caltable(self, caltable: CalTable) -> None:
        """
        A new phase table replaces the previous ones, which it was solved without
        """
        self._caltable_stack = CalStack([caltable])

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

            table = self._solve(caltable, self.solint[i])

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

            self._apply(table)

            if self.flag_dataset:
                self._flag_dataset(mode=self.flag_mode)

            self._run_imager(i)

            if self._finish_selfcal_iteration(i):
                break

            self._accept_caltable(table)
