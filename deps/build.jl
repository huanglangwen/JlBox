using BinaryProvider
using Pkg

prefix = BinaryProvider.Prefix(get(filter(!isequal("--verbose"), ARGS), 1, joinpath(@__DIR__, "usr")))
products = BinaryProvider.Product[
    BinaryProvider.FileProduct(prefix, "UManSysProp_public-1.02/umansysprop/__init__.py", :umansysprop),
]
url = "https://github.com/huanglangwen/UManSysProp_public/archive/1.02.tar.gz"
tarball_hash = "0cd723e69f770a5870ba198eb97a83617a702bb8c3dc010b16b59c1d4663698f"
unsatisfied = any(p->!BinaryProvider.satisfied(p, verbose=false), products)
if unsatisfied
    println("Can't find umansysprop in $(joinpath(@__DIR__, "usr")), downloading...")
    BinaryProvider.install(url, tarball_hash, prefix=prefix, ignore_platform=true, verbose=false)
end
#BinaryProvider.write_deps_file(joinpath(@__DIR__, "deps.jl"), products) #do not need libdl

ENV["PYTHON"]=""
Pkg.build("PyCall")

using Conda
root = Conda.ROOTENV
condaPath = Sys.iswindows() ? joinpath(root, "Scripts", "conda.exe") : joinpath(root, "bin", "conda")

run(`$(condaPath) create -n jlbox_env python=3.6 conda`)
ENV["CONDA_JL_HOME"] = "$(Conda.ROOTENV)/envs/jlbox_env"

run(`$(condaPath) install -n jlbox_env -y -c conda-forge openbabel=3.0.0`)
run(`$(condaPath) install -n jlbox_env -y Werkzeug=0.16.0 flask flask-wtf xlsxwriter`)
Pkg.build("Conda")

using Conda
println(`Current conda path is: $(Conda.ROOTENV)`)
pythonPath = Sys.iswindows() ? joinpath(Conda.PYTHONDIR, "python.exe") : joinpath(Conda.PYTHONDIR, "python")
ENV["PYTHON"]=pythonPath
Pkg.build("PyCall")
