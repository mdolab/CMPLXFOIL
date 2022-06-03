__version__ = "2.0.0"

from .pyXLIGHT import PYXLIGHT

try:
    from .postprocess import AnimateAirfoilOpt
except ImportError:
    pass
