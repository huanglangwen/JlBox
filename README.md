# JlBox

This is a rewrite version of [PyBox].

Currently only gas process implemented.

Work on Julia v0.6.

Run simulation with `include("Simulation.jl")` in julia cli or `julia -O3 Simulation.jl` in terminal.

Compared to PyBox(numba), more optimizations are (going to be) added:
- [x] constant folding for rate_values()
- [x] caching rate_values, loss_gain matrix when solving ODE
- [ ] optimize bandwidth of jacobian to use banded linear solver
- [ ] add autodiff for jacobian
- [ ] parallel version of rate_values, loss_gain and jacobian

Dependency:
- DifferentialEquations.jl
- Plots.jl
- StaticArrays.jl
Install them using `include("InstallPackages.jl")`.

[PyBox]: https://github.com/loftytopping/PyBox