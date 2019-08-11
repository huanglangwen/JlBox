using Conda
using BinaryProvider

Conda.add("openbabel",channel="conda-forge")
Conda.add("flask")
Conda.add("flask-wtf")
Conda.add("xlsxwriter")


prefix = BinaryProvider.Prefix(get(filter(!isequal("--verbose"), ARGS), 1, joinpath(@__DIR__, "usr")))
products = BinaryProvider.Product[
    BinaryProvider.LibraryProduct(prefix, "umansysprop", :umansysprop),
]
url = "https://github.com/huanglangwen/UManSysProp_public/archive/1.02.tar.gz"
tarball_hash = "0cd723e69f770a5870ba198eb97a83617a702bb8c3dc010b16b59c1d4663698f"
BinaryProvider.install(url, tarball_hash, prefix=prefix, ignore_platform=true, force=true, verbose=verbose)