using Printf
using Pkg
pkgs=["DifferentialEquations","LanguageServer","StaticArrays","CSV","DataFrames","LightXML","PyCall","Conda","Sundials","QuadGK","JLD2"]#,
for pkg in pkgs
    @printf("Installing %s\n",pkg)
    Pkg.add(pkg)
end