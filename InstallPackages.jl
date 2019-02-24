using Printf
using Pkg
pkgs=["DifferentialEquations","LanguageServer","StaticArrays","CSV","DataFrames","LightXML","PyCall","Conda"]#,
for pkg in pkgs
    @printf("Installing %s\n",pkg)
    Pkg.add(pkg)
end