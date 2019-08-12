# JlBox

This is a rewrite version of [PyBox].

Work on Julia v1.1.0.

![Structure](docs/Structure.png)

## Get Started
1. Open julia console, type `]` for Pkg management state.
2. Enter `dev https://github.com/huanglangwen/JlBox ` .
3. Enter `build JlBox` .
4. Exit julia and `cd` into the `JlBox` package folder (normally in `.julia/dev/JlBox`).
4. Run the example using `include("example/Simulation_*.jl")` in julia console.

## Features
Compared to PyBox, more optimizations are (going to be) added:
- [x] constant folding for rate_values()
- [x] caching rate_values, loss_gain matrix when solving ODE
- [x] add jacobian for gas kinetic
- [x] add jacobian for partitioning process
- [x] adjoint sensitivity analysis
- [ ] forward sensitivity analysis
- [ ] parallel version of rate_values, loss_gain and jacobian

[PyBox]: https://github.com/loftytopping/PyBox