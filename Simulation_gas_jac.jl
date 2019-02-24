include("JlBoxModule.jl")
using .Compute:run_simulation_gas_jac
using DataFrames
using CSV

@time sol,reactants2ind=run_simulation_gas_jac()
num_reactants=length(reactants2ind)
ind2reactants=Dict(reactants2ind[key]=>key for key in keys(reactants2ind))
reactants=[ind2reactants[ind] for ind in 1:num_reactants]
df=DataFrame(transpose(sol))
names!(df,[Symbol(reac) for reac in reactants])
CSV.write("data/results_jac.csv",df)
df