# JlBox

This is a rewrite version of [PyBox].

Work on Julia v1.1.0.

## Get Started
It is recommended to run JlBox in Docker:
1. Download the `Dockerfile` from this repo (in transfer branch)
2. enter the folder of `Dockerfile` in terminal and run `docker build -t jlbox .`
3. run `docker run --name=project_jlbox -it -v <some folder in your pc>:/data jlbox`
4. in the container, enter the jlbox folder: `cd /Code/Git_repos/JlBox/`
5. Run simulation with `include("Simulation_aerosol.jl")` in julia cli or `julia -O3 Simulation_aerosol.jl` in terminal.

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