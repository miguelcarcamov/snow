import os
import shutil
from pathlib import Path

from casatasks import applycal, gaincal, rmtables

from .selfcal import Selfcal


class Ampcal(Selfcal):

    def __init__(self, selfcal_object=None, input_caltable="", **kwargs):
        super(Ampcal, self).__init__(**kwargs)
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])

        self.calmode = 'a'
        self.loops = len(self.solint)
        self.imagename = self.Imager.getOutputPath()
        self.psnr_file_backup = self.Imager.output + "/psnr_amp.txt"

        if self.selfcal_object is None and self.input_caltable == "":
            print("Error, Self-cal object is Nonetype and input_caltable is an empty string")
            raise ValueError(
                "Error, Amplitude self-cal objects cannot run without an phase-cal object or input caltable"
            )
        else:
            if self.selfcal_object:
                self.input_caltable = self.selfcal_object.getCaltables()[-1]
            elif self.input_caltable != "":
                if not os.path.exists(self.input_caltable):
                    raise FileNotFoundError(
                        "The caltable " + self.input_caltable + " needs to be created"
                    )
            else:
                print("Error, Ampcal needs a non-empty list of caltables")
                raise ValueError(
                    "Error, Amplitude self-cal objects cannot run with an empty caltable list"
                )

    def run(self):
        caltable = ""
        if not self.ismodel_in_dataset():
            imagename = self.imagename + "before_apcal"
            self.Imager.run(imagename)
            print("Before amplitude self-cal: - PSNR: " + str(self.Imager.getPSNR()))
            print("Noise: " + str(self.Imager.getSTDV() * 1000.0) + " mJy/beam")
            self.write_file_backup()
            self.psnr_history.append(self.Imager.getPSNR())

            path_object = Path(self.visfile)
            current_visfile = "{0}_{2}{1}".format(
                Path.joinpath(path_object.parent, path_object.stem), path_object.suffix,
                "before_apcal"
            )
            if os.path.exists(current_visfile):
                shutil.rmtree(current_visfile)
            shutil.copytree(self.visfile, current_visfile)
        else:
            psnr_file = self.read_first_line_file_backup()
            self.psnr_history.append(psnr_file)
            self.delete_last_lines()

        for i in range(0, self.loops):
            caltable = self.output_caltables + 'ampcal_' + str(i)
            self.caltables.append(caltable)
            rmtables(caltable)

            self.set_attributes_from_dicts(i)

            gaincal(
                vis=self.visfile,
                field=self.Imager.getField(),
                caltable=caltable,
                spw=self.Imager.getSpw(),
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
            self.save_selfcal(caltable_version=versionname, overwrite=True)
            self.caltables_versions.append(versionname)
            applycal(
                vis=self.visfile,
                spw=self.Imager.getSpw(),
                spwmap=[self.spwmap, self.spwmap],
                field=self.Imager.getField(),
                gaintable=[self.input_caltable, caltable],
                gainfield='',
                calwt=False,
                flagbackup=False,
                interp=self.interp,
                applymode=self.applymode
            )

            if self.flag_dataset_bool:
                self.flag_dataset(mode=self.flag_mode)

            imagename = self.imagename + '_a' + str(i)

            self.Imager.run(imagename)

            self.psnr_history.append(self.Imager.getPSNR())

            print("Solint: " + str(self.solint[i]) + " - PSNR: " + str(self.psnr_history[i]))
            print("Noise: " + str(self.Imager.getSTDV() * 1000.0) + " mJy/beam")

            path_object = Path(self.visfile)
            current_visfile = "{0}_{2}{1}".format(
                Path.joinpath(path_object.parent, path_object.stem), path_object.suffix,
                "amp" + str(i)
            )
            if self.restore_PSNR:
                if i > 0:
                    if self.psnr_history[-1] <= self.psnr_history[-2]:
                        self.restore_selfcal(caltable_version=self.caltables_versions[-1])
                        self.psnr_history.pop()
                        self.caltables.pop()
                        print(
                            "PSNR decreasing in this solution interval - restoring to last MS and exiting loop"
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

                elif self.selfcal_object:
                    if self.psnr_history[i] <= self.selfcal_object.getPSNRHistory()[-1]:
                        self.restore_selfcal(
                            caltable_version=self.selfcal_object.getCaltablesVersions()[-1]
                        )
                        self.psnr_history.pop()
                        self.caltables.pop()
                        self.caltables = self.selfcal_object.getCaltables()
                        self.psnr_history = self.selfcal_object.getPSNRHistory()
                        self.caltables_versions = self.selfcal_object.getCaltablesVersions()
                        print(
                            "PSNR decreasing or equal in this solution interval - restoring to last MS and "
                            "exiting loop"
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
