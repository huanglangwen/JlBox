# JlBox

This is a rewrite version of [PyBox].

Work on Julia v1.1.0.

## Get Started
WIP, see `transter` branch.

## Features
Compared to PyBox, more optimizations are (going to be) added:
- [x] constant folding for rate_values()
- [x] caching rate_values, loss_gain matrix when solving ODE
- [x] add jacobian for gas kinetic
- [x] add jacobian for partitioning process
- [x] adjoint sensitivity analysis
- [ ] forward sensitivity analysis
- [ ] parallel version of rate_values, loss_gain and jacobian

## Dependency
- DifferentialEquations.jl
- StaticArrays.jl
- DataFrames.jl
- CSV.jl
- Conda.jl
- PyCall.jl
- LightXML.jl
- QuadGK.jl
Install them using `include("InstallPackages.jl")`.

[PyBox]: https://github.com/loftytopping/PyBox