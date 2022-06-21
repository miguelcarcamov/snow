from .imager import Imager


class WSClean(Imager):

    def __init__(self, **kwargs):
        super(WSClean, self).__init__(**kwargs)
        initlocals = locals()
        initlocals.pop('self')
        for a_attribute in initlocals.keys():
            setattr(self, a_attribute, initlocals[a_attribute])
        self.__dict__.update(kwargs)

    def run(self, imagename=""):
        return
