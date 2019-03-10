using Printf
using Pkg
pkgs=["DifferentialEquations","StaticArrays"]#,
for pkg in pkgs
    @printf("Installing %s\n",pkg)
    Pkg.add(pkg)
end