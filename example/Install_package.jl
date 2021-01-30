using Pkg
pkgs = ["DataFrames", "CSV", "ArgParse", "Sundials", "OrdinaryDiffEq"]
for pkg in pkgs
    Pkg.add(pkg)
end