"""Draw the figures of the manual's Experimental geometry page.

    PYTHONPATH=src python3.14 tools/geometry_figures.py [folder]

writes, into docs/assets/ (or ``folder``):

- geometry-ibm.png: GEOMETRY IBM as in a chamber -- the sample turning about a
  vertical tilt axis, with beam, detector and normal in the red plane: for
  THETA +25 and -25, the 3D view, the top view and the RUMP commands.
- geometry-cornell.png: GEOMETRY CORNELL as in a chamber -- the sample
  turning about a vertical tilt axis like a door, with the plane of beam and
  normal (red), the plane of beam and detector (orange), the top and front
  views and the RUMP commands.
- geometry-general.png: GEOMETRY GENERAL in the same perspective -- the
  normal rising out of the red plane, so PSI takes a value no rule computes.

Each figure is its own module under tools/geometry/, runnable on its own
(``tools/geometry/cornell.py out.png``); tools/geometry/common.py holds the
angles they are drawn at (THETA 25, PHI 35), the colour code and the arrows
and arcs. The exit angles in the command boxes are pyRUMP's own
(``pyrump.model.geometry.Geometry``), so the figures can't drift from the
code. They say "detector", not "RBS detector": the same angles describe the
PIXE detector.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent / "geometry"))

import cornell  # noqa: E402
import general  # noqa: E402
import ibm  # noqa: E402

FIGURES = {
    "geometry-ibm.png": ibm.figure,
    "geometry-cornell.png": cornell.figure,
    "geometry-general.png": general.figure,
}


def main(folder: Path) -> None:
    for name, draw in FIGURES.items():
        draw(folder / name)
        print(f"wrote {folder / name}")


if __name__ == "__main__":
    default = Path(__file__).resolve().parents[1] / "docs" / "assets"
    main(Path(sys.argv[1]) if len(sys.argv) > 1 else default)
