__version__ = "2.1.1"

from .CMPLXFOIL import CMPLXFOIL

try:
    from .postprocess import AnimateAirfoilOpt
except ImportError:
    pass
