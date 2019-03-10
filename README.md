# JlBox

This is a minimal version of JlBox for Reaction Network mechanisms.

Work on Julia v1.1.0.

## Get Started
1. run `InstallPackages.jl` at the first time
2. Change `Configure_gas.jl` if you like
3. run `Simulation_gas.jl`


## Features
Compared to PyBox, more optimizations are (going to be) added:
- [x] constant folding for rate_values()
- [x] caching rate_values, loss_gain matrix when solving ODE
- [x] add jacobian for gas kinetic
- [ ] parallel version of rate_values, loss_gain and jacobian

## Dependency
- DifferentialEquations.jl
- StaticArrays.jl
Install them using `include("InstallPackages.jl")`.

[PyBox]: https://github.com/loftytopping/PyBox