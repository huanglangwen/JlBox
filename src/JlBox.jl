module JlBox
using FiniteDiff
using ForwardDiff
using ForwardDiff:JacobianConfig
using DiffEqBase
using OrdinaryDiffEq
using StaticArrays
using SparseArrays
using QuadGK
using Serialization
using PyCall
using LightXML
using LinearAlgebra
using DiffEqCallbacks
using DiffEqOperators
using Logging
using Requires

include(joinpath(@__DIR__,"../deps/deps.jl"))

function __init__()
    check_deps()
    @require Plots="91a5bcdd-55d7-5caf-9e0b-520d859cae80" include("Plot.jl")
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