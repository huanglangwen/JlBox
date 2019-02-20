#using Printf
pkgs=["DifferentialEquations","LanguageServer","StaticArrays","CSV","DataFrames"]#,
for pkg in pkgs
    @printf("Installing %s\n",pkg)
    Pkg.add(pkg)
end