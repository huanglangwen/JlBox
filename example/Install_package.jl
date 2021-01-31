using Pkg
pkgs = ["DataFrames", "CSV", "ArgParse", "Sundials", "OrdinaryDiffEq", "Plots"]
for pkg in pkgs
    Pkg.add(pkg)
end