using Pkg
pkgs = ["Sundials", "DataFrames", "CSV", "ArgParse", "OrdinaryDiffEq"]
for pkg in pkgs
    Pkg.add(pkg)
end