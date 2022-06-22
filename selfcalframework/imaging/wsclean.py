from .imager import Imager


class WSClean(Imager):

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def run(self, imagename=""):
        pass
