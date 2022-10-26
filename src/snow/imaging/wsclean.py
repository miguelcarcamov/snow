from dataclasses import dataclass
from .imager import Imager


@dataclass(init=True, repr=True)
class WSClean(Imager):

    def run(self, imagename=""):
        """
        Method that runs WSClean

        Parameters
        ----------
        imagename :
            The absolute path to the output image name
        """
        raise NotImplementedError("This function has not been implemented yet")
