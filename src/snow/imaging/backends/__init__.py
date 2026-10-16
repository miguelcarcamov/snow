"""
Built-in imagers that follow the contract of snow.imaging.contract (#23). Optional backends, such as
Pyralysis (#24), are imported lazily so their dependencies stay optional.
"""
from .tclean import Tclean
