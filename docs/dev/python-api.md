## Python API

The CLI is a thin wrapper; the library is the primary interface.

```python
from pathlib import Path
import pyrump
from pyrump.atomic.density import DensityTable
from pyrump.atomic.tables import PeriodicTable
from pyrump.io.kalbitzer import parse_kalbitzer
from pyrump.model.detector import Measurement
from pyrump.model.geometry import Geometry
from pyrump.model.spectrum import Calibration
from pyrump.sim.engine import Beam, UniformSample, simulate
from pyrump.stopping.kalbitzer import KalbitzerStopping
from pyrump.stopping.registry import StoppingRegistry
from pyrump.stopping.ziegler import ZieglerStopping

DATA = Path(pyrump.__file__).parent / "data"   # bundled with the package

table = PeriodicTable.load(DATA / "atom4.dat", DATA / "pscoef.dat")
registry = StoppingRegistry(
    table.elements,
    kalbitzer=KalbitzerStopping(parse_kalbitzer(f"{DATA}/newstop.kal"), table.elements),
    ziegler=ZieglerStopping(table.elements),
)

spectrum = simulate(
    UniformSample(
        thicknesses=[1000.0],        # 1e15 atoms/cm^2
        element_z=[14],              # silicon
        compositions=[[1.0]],
    ),
    Beam(e0_MeV=2.0, z=2, mass=4.0026),
    Geometry(theta=0.0, phi=10.0),   # phi = 180 - scattering angle
    registry,
    table,
    Calibration(kevch=5.0, kev0=0.0, npt=1024),
    Measurement(omega_msr=1.0, charge_uC=10.0, fwhm_keV=15.0),
)

print(spectrum.total(), "counts")
```

Building the registry takes a moment; **build it once and reuse it**, especially
when fitting.

### Reading and writing files

```python
from pyrump.io.rbs import read_rbs, write_rbs
from pyrump.io.ascii import read_ascii, write_ascii

measured = read_rbs("data.rbs")
measured.counts          # np.ndarray
measured.calibration     # keV/channel, offset
measured.geometry        # angles, geometry convention
measured.e0_MeV, measured.zbeam, measured.mbeam
```

