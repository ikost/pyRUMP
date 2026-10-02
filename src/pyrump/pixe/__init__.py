"""PIXE: particle-induced X-ray emission, simulated from the same sample as RBS.

See ``docs/manual/pixe.md`` for the design.
"""

from .data import PixeData
from .detector import DEFAULT_CALIBRATION, Absorber, PixeDetector

__all__ = ["DEFAULT_CALIBRATION", "Absorber", "PixeData", "PixeDetector"]
