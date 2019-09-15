using Conda
using BinaryProvider
using Pkg

if Conda.exists("openbabel==2.4.1")#Conda.exists("openbabel") doesn't work now!
    Conda.add("openbabel",channel="conda-forge")
end
Conda.add("flask")
Conda.add("flask-wtf")
Conda.add("xlsxwriter")

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
