## Things that will catch you out

**`phi` is not the scattering angle.** It is 180° minus it. A detector at 170°
means `phi = 10`. Use `geometry.scattering_angle` for the physical value.

**Fitting windows must cover channels where the model has counts.** Poisson
likelihood is undefined where the model predicts zero, so those channels
contribute neither to χ² nor to the gradient. A window reaching past the
spectrum silently throws away most of its own evidence, and a parameter can sit
motionless while the fit reports success:

```
warning: 81 of 81 windowed channels had zero predicted counts
  thickness[0]                     995.54     ← did not move
```

pyRUMP reports the count and warns; RUMP's manual mentions it in one line.
Narrow the window and refit.

**Thickness is areal density**, 10¹⁵ atoms/cm². Converting to nanometres needs
an assumed atomic density.

**Straggling is off by default**, matching RUMP. Set `straggle=1.0` for the
Bohr value. Note the in/out combination is an approximation and there is no Chu
correction.

**Absolute yields depend on charge, solid angle and efficiency**, which are
rarely known to better than a few percent. Use a normalisation window rather
than trusting them — the worked example above is 9× out for exactly this reason.

**Neighbouring heavy elements are not separable.** In and Sn differ by 0.6
channels at 3 MeV; fit their ratio, not each independently.

**Build the stopping registry once.** It refits polynomials per beam energy, so
recreating it inside a fit loop is slow for no reason.

