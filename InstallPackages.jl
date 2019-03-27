using Printf
using Pkg
pkgs=["DifferentialEquations","LanguageServer","StaticArrays","CSV","DataFrames","LightXML","PyCall","Conda","Sundials","QuadGK","BSON"]#,
for pkg in pkgs
    @printf("Installing %s\n",pkg)
    Pkg.add(pkg)
end