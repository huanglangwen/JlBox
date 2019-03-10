using Printf
using Pkg
pkgs=["DifferentialEquations","StaticArrays","TimerOutputs","OrdinaryDiffEq"]#,
for pkg in pkgs
    @printf("Installing %s\n",pkg)
    Pkg.add(pkg)
end