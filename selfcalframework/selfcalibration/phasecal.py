import os
import shutil
from pathlib import Path

from casatasks import applycal, gaincal, rmtables

from .selfcal import Selfcal


class Phasecal(Selfcal):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

        self._calmode = 'p'
        self._loops = len(self.solint)

    def run(self):
        caltable = "before_selfcal"
        self._save_selfcal(caltable_version=caltable, overwrite=True)
        self.caltables_versions.append(caltable)
        if not self._ismodel_in_dataset() and self.previous_selfcal is None:
            imagename = self._image_name + '_original'
            self.imager.run(imagename)
            print("Original: - PSNR: " + str(self.imager.psnr))
            print("Noise: " + str(self.imager.stdv * 1000.0) + " mJy/beam")
            self._psnr_history.append(self.Imager.getPSNR())
        elif self.previous_selfcal is not None:
            self.psnr_history.append(self.previous_selfcal._psnr_history[-1])

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
                calmode=self.calmode,
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

            imagename = self.imagename + '_ph' + str(i)

            self.imager.run(imagename)

            self.psnr_history.append(self.imager.psnr)

            print("Solint: " + str(self.solint[i]) + " - PSNR: " + str(self._psnr_history[-1]))
            print("Noise: " + str(self.imager.stdv * 1000.0) + " mJy/beam")
            path_object = Path(self.visfile)
            current_visfile = "{0}_{2}{1}".format(
                Path.joinpath(path_object.parent, path_object.stem), path_object.suffix,
                "ph" + str(i)
            )

            if os.path.exists(current_visfile):
                shutil.rmtree(current_visfile)
            shutil.copytree(self.visfile, current_visfile)

            if self.restore_psnr:
                if len(self._psnr_history) > 1:
                    if self._psnr_history[-1] <= self._psnr_history[-2]:
                        self._restore_selfcal(caltable_version=self._caltables_versions[i])
                        self._psnr_history.pop()
                        self._caltables.pop()
                        print(
                            "PSNR decreasing or equal in this solution interval - restoring to last MS and exiting loop"
                        )
                        break
                    else:
                        print(
                            "PSNR improved on iteration {0} - Copying measurement set files...".
                            format(i)
                        )
                        if os.path.exists(current_visfile):
                            shutil.rmtree(current_visfile)
                        shutil.copytree(self.visfile, current_visfile)
                        self.visfile = current_visfile
                        self.imager.inputvis = current_visfile
