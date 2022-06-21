import os
import shutil
from pathlib import Path

from casatasks import applycal, gaincal, rmtables

from .selfcal import Selfcal


class Phasecal(Selfcal):

    def __init__(self, input_caltable="", **kwargs):
        super(Phasecal, self).__init__(**kwargs)
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])

        self.calmode = 'p'
        self.loops = len(self.solint)
        self.imagename = self.Imager.getOutputPath()
        self.psnr_file_backup = self.Imager.output + "/psnr_ph.txt"

    def run(self):
        caltable = "before_selfcal"
        self.save_selfcal(caltable_version=caltable, overwrite=True)
        self.caltables_versions.append(caltable)
        if not self.ismodel_in_dataset():
            imagename = self.imagename + '_original'
            self.Imager.run(imagename)
            print("Original: - PSNR: " + str(self.Imager.getPSNR()))
            print("Noise: " + str(self.Imager.getSTDV() * 1000.0) + " mJy/beam")
            self.write_file_backup()
            self.psnr_history.append(self.Imager.getPSNR())
        else:
            psnr_file = self.read_first_line_file_backup()
            self.psnr_history.append(psnr_file)
            self.delete_last_lines()

            path_object = Path(self.visfile)
            current_visfile = "{0}_{2}{1}".format(
                Path.joinpath(path_object.parent, path_object.stem), path_object.suffix, "original"
            )
            if os.path.exists(current_visfile):
                shutil.rmtree(current_visfile)
            shutil.copytree(self.visfile, current_visfile)

        for i in range(0, self.loops):
            caltable = 'pcal' + str(i)
            self.caltables.append(caltable)
            rmtables(caltable)

            self.set_attributes_from_dicts(i)

            gaincal(
                vis=self.visfile,
                caltable=caltable,
                field=self.Imager.getField(),
                spw=self.Imager.getSpw(),
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

            self.plot_selfcal(
                caltable,
                xaxis="time",
                yaxis="phase",
                iteration="antenna",
                subplot=[4, 2],
                plotrange=[0, 0, -180, 180],
                want_plot=self.want_plot
            )

            versionname = 'before_phasecal_' + str(i)
            self.save_selfcal(caltable_version=versionname, overwrite=True)
            self.caltables_versions.append(versionname)

            applycal(
                vis=self.visfile,
                field=self.Imager.getField(),
                spw=self.Imager.getSpw(),
                spwmap=self.spwmap,
                gaintable=[caltable],
                gainfield='',
                calwt=False,
                flagbackup=False,
                interp=self.interp,
                applymode=self.applymode
            )

            if self.flag_dataset_bool:
                self.flag_dataset(mode=self.flag_mode)

            imagename = self.imagename + '_ph' + str(i)

            self.Imager.run(imagename)

            self.psnr_history.append(self.Imager.getPSNR())

            print("Solint: " + str(self.solint[i]) + " - PSNR: " + str(self.psnr_history[-1]))
            print("Noise: " + str(self.Imager.getSTDV() * 1000.0) + " mJy/beam")
            path_object = Path(self.visfile)
            current_visfile = "{0}_{2}{1}".format(
                Path.joinpath(path_object.parent, path_object.stem), path_object.suffix,
                "ph" + str(i)
            )
            if os.path.exists(current_visfile):
                shutil.rmtree(current_visfile)
            shutil.copytree(self.visfile, current_visfile)
            if self.restore_PSNR:
                if len(self.psnr_history) > 1:
                    if self.psnr_history[-1] <= self.psnr_history[-2]:
                        self.restore_selfcal(caltable_version=self.caltables_versions[i])
                        self.psnr_history.pop()
                        self.caltables.pop()
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
                        self.Imager.inputvis = current_visfile
                        self.write_file_backup()
