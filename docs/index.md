# pyRUMP

<div class="term-window">
  <div class="term-titlebar">
    <span class="term-dot term-dot--red"></span>
    <span class="term-dot term-dot--yellow"></span>
    <span class="term-dot term-dot--green"></span>
  </div>
  <div class="term-body">
    <pre class="ascii-banner">
########  ##    ## ########  ##     ## ##     ## ########  
##     ##  ##  ##  ##     ## ##     ## ###   ### ##     ## 
##     ##   ####   ##     ## ##     ## #### #### ##     ## 
########     ##    ########  ##     ## ## ### ## ########  
##           ##    ##   ##   ##     ## ##     ## ##        
##           ##    ##    ##  ##     ## ##     ## ##        
##           ##    ##     ##  #######  ##     ## ##        
    </pre>
    <p class="term-prompt"><span class="term-caret">$</span> pyrump</p>
  </div>
</div>

An AI-powered Python reimplementation of **RUMP**, the Rutherford backscattering
spectrometry (RBS) simulation and analysis package originally written by
L. R. Doolittle and M. O. Thompson at Cornell. pyRUMP installs via a standard
`pip install` and runs on Windows, macOS, and Linux.

The original RUMP is ~22k lines of C code from the late 1980s (archived for historical
reference at [github.com/ikost/Rump-and-Genplot-Archive](https://github.com/ikost/Rump-and-Genplot-Archive)), with an
HTML manual and no active support; the distribution was previously available via
`genplot.com`, now abandoned and no longer resolving.
pyRUMP reproduces the simulation physics and data-analysis flow as a tested, importable
Python library, with both a classical RUMP CLI and a Python API. Plotting uses matplotlib;
fitting uses scipy/numpy.

Each of these was a deliberate choice. The original C code needed a separate makefile per
platform (`makeaix`, `makelnx`, `makeosx`, `makesgi`, `makesolaris`, ...); Python collapses
that into one `pip install` that runs unmodified on Windows, macOS, and Linux. scipy and
numpy provide compiled, vectorized numerical routines, so the least-squares fitting and
data analysis stay fast without any hand-written C. And matplotlib's backend system
replaces the dozen-plus hardware-specific printer and plotter drivers GENPLOT used to ship
(`deskjet`, `epson`, `hpgl`, `laserjet`, `postdrv`, `tektrx`, `xdriver`, and more) with a
single plotting API that targets any output device.

The original RUMP ecosystem included **RUMPX**, an X-Window GUI popular among Windows
users two decades ago. pyRUMP doesn't have a GUI counterpart yet — the author prefers
CLI tools — but a Python equivalent could be added later if there's real demand for it.

![Measured spectrum compared against a pyRUMP simulation, showing a Si / Ru / Mn2.74Pt1 / Ru stack](assets/mnpt-fit.png)
*Measured data (black) vs. a pyRUMP simulation (red) of a Si / Ru / Mn<sub>2.74</sub>Pt<sub>1</sub> / Ru stack.*

## Where to go next

- **[Getting started](getting-started/index.md)** — install pyRUMP and run your first simulation
- **[Manual](manual/physics.md)** — what the forward model computes, the interactive shell, CLI, and sample-description reference
- **[Development](dev/python-api.md)** — the Python API, contributing, validation against the original, and porting notes
