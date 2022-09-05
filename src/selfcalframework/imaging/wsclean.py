from .imager import Imager


class WSClean(Imager):

    def __init__(self, **kwargs):
        """
        WSClean object

        Parameters
        ----------
        kwargs :
            General imager arguments
        """
        super().__init__(**kwargs)

    def run(self, imagename=""):
        """
        Method that runs WSClean

        Parameters
        ----------
        imagename :
            The absolute path to the output image name
        """
        raise NotImplementedError("This function has not been implemented yet")
