module JlBox
using FiniteDiff
using ForwardDiff
using ForwardDiff:JacobianConfig
using OrdinaryDiffEq
using StaticArrays
using SparseArrays
using Printf
using QuadGK
using Serialization
using PyCall
using LightXML
using LinearAlgebra

include("../deps/deps.jl")

function __init__()
    check_deps()
end

include("Configure.jl")
include("Optimize.jl")
include("Parse_eqn.jl")
include("Size_Distributions.jl")
include("PropertyCalculation.jl")
include("Partitioning.jl")
include("Jacobian.jl")
include("Sensitivity.jl")
include("Prepare.jl")
include("Compute.jl")

export run_simulation_gas, run_simulation_aerosol, run_simulation_aerosol_adjoint
end