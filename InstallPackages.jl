using Printf
using Pkg
pkgs=["DifferentialEquations","StaticArrays","TimerOutputs"]#,
for pkg in pkgs
    @printf("Installing %s\n",pkg)
    Pkg.add(pkg)
end