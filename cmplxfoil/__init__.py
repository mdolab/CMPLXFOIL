__version__ = "2.0.4"

from .CMPLXFOIL import CMPLXFOIL

try:
    from .postprocess import AnimateAirfoilOpt
except ImportError:
    pass
