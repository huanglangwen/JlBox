# JlBox

This is a rewrite version of [PyBox].

Work on Julia v0.6.

Run simulation with `include("Simulation_aerosol.jl")` in julia cli or `julia -O3 Simulation_aerosol.jl` in terminal.

Compared to PyBox(numba), more optimizations are (going to be) added:
- [x] constant folding for rate_values()
- [x] caching rate_values, loss_gain matrix when solving ODE
- [ ] add jacobian for gas kinetic and partitioning process
- [ ] parallel version of rate_values, loss_gain and jacobian

Dependency:
- DifferentialEquations.jl
- StaticArrays.jl
- DataFrames.jl
- CSV.jl
Install them using `include("InstallPackages.jl")`.

[PyBox]: https://github.com/loftytopping/PyBox