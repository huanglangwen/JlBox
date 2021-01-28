using Pkg
pkgs = ["DataFrames", "CSV", "ArgParse"]
for pkg in pkgs
    Pkg.add(pkg)
end