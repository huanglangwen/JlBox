# JlBox

The `JlBox` is a julia package that simulates the evolution of chemicals in the atmosphere using
box model where advection effect is ignored. It is heavily inspired by Dr. David Topping's [PyBox]
and the two models could produce identical results.

This package works on Julia v1.1 .

## Get Started

### Option 1: Running on an existing julia environment
This should work on Windows, MacOS and Linux, but I only tested the Linux one.
1. Open julia console, type `]` for Pkg management state.
2. Enter `dev https://github.com/huanglangwen/JlBox ` .
3. Enter `build JlBox` .
4. Exit julia and `cd` into the `JlBox` package folder (normally in `.julia/dev/JlBox`).
4. Run the example using `include("example/Simulation_*.jl")` in julia console.

### Option 2: Running on Docker
This method should guarantee to run since it creates a standard Linux environment.
1. Download the `Dockerfile` from this repo (in transfer branch)
2. enter the folder of `Dockerfile` in terminal and run `docker build -t jlbox .`
3. run `docker run --name=project_jlbox -it jlbox`
4. in the container, enter the jlbox folder: `cd /root/.julia/dev/JlBox`
5. Run simulation with `include("example/Simulation_*.jl")` in julia cli.

### Custom Simulation
The `JlBox` package exposes functions like `run_simulation_*(config)`. Indicated by the function
name, there are three simulation types: gas only simulation (only gas kinetics process), gas-aerosol simulation (gas kinetics+gas-aerosol partitioning) and adjoint sensitivity analysis for the gas-aerosol
simulation. Users could tweak some parameters by changing the `config` as is illustrated in the examples.

## Features
Compared to PyBox, more optimizations are (going to be) added:
- [x] constant folding for rate_values()
- [x] caching rate_values, loss_gain matrix when solving ODE
- [x] add jacobian for gas kinetic
- [x] add jacobian for partitioning process
- [x] adjoint sensitivity analysis
- [ ] sparse auto-differentiation
- [ ] native ode solvers
- [ ] forward sensitivity analysis
- [ ] parallel version of rate_values, loss_gain and jacobian

## Internals
![Structure](docs/Structure.png)

[PyBox]: https://github.com/loftytopping/PyBox