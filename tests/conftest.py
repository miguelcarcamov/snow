"""Shared test configuration.

CASA is not installable everywhere (wheels exist only for some platforms and Python versions), so
CASA-free unit tests stub ``casatasks``/``casatools`` when they are missing. Tests that need real
CASA must be marked ``@pytest.mark.casa`` and live under ``tests/integration``.
"""
import importlib.util
import sys
from unittest import mock

for _module in ("casatasks", "casatools"):
    if importlib.util.find_spec(_module) is None:
        sys.modules[_module] = mock.MagicMock(name=_module)
