#using Printf
pkgs=["Plots","DifferentialEquations","LanguageServer","StaticArrays"]#,
for pkg in pkgs
    @printf("Installing %s\n",pkg)
    Pkg.add(pkg)
end