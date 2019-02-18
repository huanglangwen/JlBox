include("JlBoxModule.jl")
using Compute:run_simulation_aerosol
using DataFrames
using CSV

@time sol,reactants2ind=run_simulation_aerosol()
size(sol)
