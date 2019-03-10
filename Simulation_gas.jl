include("JlBoxModule.jl")
using .Compute:run_simulation_gas
#using DataFrames
#using CSV

@time sol,reactants2ind=run_simulation_gas()

#using Plots
#plot(log10.(sol[reactants2ind["BUT1ENE"],:]))
#plot!(log10.(sol[reactants2ind["C2H5O2"],:]))