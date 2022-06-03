__version__ = "2.0.0"

from .CMPLXFOIL import CMPLXFOIL

try:
    from .postprocess import AnimateAirfoilOpt
except ImportError:
    pass
