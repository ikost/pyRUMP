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
The original RUMP ecosystem included **RUMPX**, an X-Window GUI popular among Windows
users two decades ago. pyRUMP doesn't have a GUI counterpart yet — the author prefers
CLI tools — but a Python equivalent could be added later if there's real demand for it.

![Measured spectrum compared against a pyRUMP simulation, showing a Pt/MnGeN/Pt stack on Si](assets/mnpt-fit.png)
*Measured data (black) vs. a pyRUMP simulation (red) of a Si / Ru / Mn<sub>2.74</sub>Pt<sub>1</sub> / Ru stack.*

## Where to go next

- **[Getting started](getting-started/index.md)** — install pyRUMP and run your first simulation
- **[Physics and simulation](physics/index.md)** — what the forward model computes, and how well it agrees with the original
- **[Manual](manual/shell.md)** — the interactive shell, CLI, and sample-description reference
- **[Development](dev/python-api.md)** — the Python API, contributing, and porting notes
