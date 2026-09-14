## Worked examples

Every number below was produced by running the code, not written from memory.
For what the simulation actually computes, see
[How the simulation works](../physics/index.md).

### Identifying what is in a sample

Given an unknown spectrum, the first question is which elements are present.
Each element's **surface edge** sits at *K*·*E*₀, so predicted edge positions
identify the peaks.

Using `2A.rbs`, one of the files shipped with RUMP (see
[Licensing and provenance](../dev/about.md#licensing-and-provenance) — it isn't redistributed
with pyRUMP, so point `PYRUMP_DATA` at your own copy of the legacy `rump/data/`
tree to reproduce this):

```python
import numpy as np
from pyrump.io.rbs import read_rbs
from pyrump.physics.kinematics import kinematic_factor

s = read_rbs("C-code/rump/data/Fixed/2A.rbs")
print(s.identifier)
print(f"{s.e0_MeV} MeV, Z={s.zbeam}, scattering angle {s.geometry.scattering_angle}")
```

```
Binghampton_target_02A.RBS  RBS LT =  905.98 RT  962.42
3.0 MeV, Z=1, scattering angle 160.0
```

So: 3 MeV protons at 160°, 7.815 keV/channel with a 65.6 keV offset. Now
predict where each candidate element's edge would fall:

```python
E0 = s.e0_MeV * 1000
for symbol in ("C", "O", "Si", "Ti", "Fe", "In", "Sn", "Au"):
    element = table.by_symbol(symbol)
    mass = max(element.isotopes, key=lambda i: i.fraction).mass
    K = kinematic_factor(s.mbeam, mass, s.geometry.scattering_angle)
    print(f"{symbol:3s} K={K:.4f}  E={K*E0:7.1f} keV  channel {s.calibration.channel_of(K*E0):6.1f}")
```

```
C   K=0.7213  E= 2164.0 keV  channel  268.5
O   K=0.7829  E= 2348.6 keV  channel  292.1
Si  K=0.8695  E= 2608.5 keV  channel  325.4
Ti  K=0.9217  E= 2765.0 keV  channel  345.4
Fe  K=0.9325  E= 2797.4 keV  channel  349.6
In  K=0.9665  E= 2899.6 keV  channel  362.6
Sn  K=0.9679  E= 2903.7 keV  channel  363.2
Au  K=0.9803  E= 2941.0 keV  channel  367.9
```

The measured spectrum has falling edges at channels **267, 292 and 363**, which
match **carbon, oxygen, and indium/tin**. Indium and tin are 0.6 channels apart
here and cannot be separated — a general limitation for neighbouring heavy
elements, and the reason a fit constrains their *ratio* rather than resolving
them independently.

That composition — In, Sn, O over a C/O/H substrate — is indium tin oxide on a
polymer, which is exactly what `ITO.lcm` in the same directory describes.

### Simulating a known structure

RUMP ships both the measurement and a matching sample description, so we can
simulate one against the other:

```python
from pyrump.script.lcm import read_lcm, to_sample
from pyrump.sim.engine import Beam, simulate

densities = DensityTable.load(f"{DATA}/density.tab")
observed  = read_rbs(f"{DATA}/Fixed/2A.rbs")
sample    = to_sample(read_lcm(f"{DATA}/Fixed/ITO.lcm"), table, densities)

simulated = simulate(
    sample,
    Beam(e0_MeV=observed.e0_MeV, z=observed.zbeam, mass=observed.mbeam),
    observed.geometry, registry, table,
    observed.calibration, observed.measurement,
)
```

Comparing yields in each edge region:

| Region | Channels | Measured | Simulated |
|---|---|---|---|
| C edge | 255–270 | 266 902 | 27 205 |
| O edge | 280–295 | 63 160 | 8 150 |
| In/Sn edge | 350–366 | 35 261 | 2 273 |
| **total** | | **3 489 801** | **380 110** |

The **structure is right** — the edges land in the right channels and the
relative intensities are close — but the absolute yield is **9.2× low**.

That is not a simulation error: running the same case through the original C
gives 381 120 counts, agreeing with pyRUMP to 2.6e-3. Both codes say the same
thing, so the discrepancy lives in the measurement's normalisation — the
actual collected charge, solid angle, or detector efficiency differs from the
values recorded in the file. **This is the normal situation in RBS**, which is
why RUMP has both a `CORR` factor and a normalisation window. Rather than
trusting the charge integration, you fit the scale:

```python
from pyrump.fit.windows import Window, WindowSet

windows = WindowSet(
    error=[Window(255, 370)],          # fit over the interesting region
    normalisation=Window(255, 370),    # and let the scale float
)
```

The normalisation window forces the total counts over that range to agree by
scaling the data, before χ² is evaluated — so a charge-integration error stops
biasing the fitted thicknesses.

### Fitting a thickness

Simulate a 2400 Å silicon layer, then recover it from a 2000 Å starting guess.

```bash
pyrump simulate truth.lcm  --energy 2.0 -o data.rbs
pyrump fit      start.lcm data.rbs --energy 2.0 --vary thickness:0 --window 190 226
```

```
reduced chi-square 0.0000 on 36 dof
10 evaluations, `xtol` termination condition is satisfied.
  thickness[0]                     1194.6  +/- 2.7494
```

2400 Å of silicon is 2400 × 0.4977 = **1194.5** in 10¹⁵ atoms/cm², so the fit
recovers it to better than 0.1%. Note that thickness is reported in **areal
density**, not Ångström — RBS measures atoms per unit area, and converting to
a physical thickness needs an assumed density, a separate and often less
certain quantity.

The same fit from Python, with two parameters:

```python
from pyrump.fit.lm import fit
from pyrump.fit.parameters import FitInputs, thickness, parameter
from pyrump.fit.windows import Window, WindowSet

inputs = FitInputs(sample=sample, beam=beam, geometry=geometry,
                   calibration=calibration, measurement=measurement)

result = fit(
    lambda i: simulate(i.sample, i.beam, i.geometry, registry, table,
                       i.calibration, i.measurement).counts,
    measured_counts, inputs,
    [thickness(0), parameter("fwhm")],
    windows=WindowSet(error=[Window(200, 260)]),
)

print(result.parameters)      # {'thickness[0]': 1197.0, 'fwhm': 14.77}
print(result.uncertainties)
print(result.correlation)     # parameters are rarely independent
```

Always look at the correlation matrix. Thickness and resolution were −0.35
correlated here; strongly correlated parameters mean the data does not
constrain them separately, however tight the individual error bars look.

