module JlBox
using FiniteDiff
using ForwardDiff
using ForwardDiff:JacobianConfig
using DiffEqBase
using OrdinaryDiffEq
using StaticArrays
using SparseArrays
using Printf
using QuadGK
using Serialization
using PyCall
using LightXML
using LinearAlgebra
using DiffEqCallbacks
using DiffEqOperators
using Logging

include("../deps/deps.jl")

function __init__()
    check_deps()
end

include("Constant.jl")
include("Configure.jl")
include("Optimize.jl")
include("Parse_eqn.jl")
include("Size_Distributions.jl")
include("PropertyCalculation.jl")
include("Partitioning.jl")
include("Jacobian.jl")
include("Sensitivity.jl")
include("Prepare.jl")
include("RHS.jl")
include("Preconditioning.jl")
include("Simulation.jl")
include("PostProcessing.jl")

export run_simulation, run_simulation_aerosol_adjoint
end