const umansysprop = joinpath(dirname(@__FILE__), "usr/UManSysProp_public-1.02/umansysprop/__init__.py")
function check_deps()
    global umansysprop
    if !isfile(umansysprop)
        error("$(umansysprop) does not exist, Please re-run Pkg.build(\"JlBox\"), and restart Julia.")
    end

end
