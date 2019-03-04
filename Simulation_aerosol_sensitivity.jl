include("JlBoxModule.jl")
include("Configure_aerosol.jl")
using .Compute:run_simulation_aerosol_sensitivity
using DataFrames
using CSV

@time lambda_sol,reactants2ind,num_reactants=run_simulation_aerosol_sensitivity(linsolver=:Dense)

ind2reactants=Dict(reactants2ind[key]=>key for key in keys(reactants2ind))
reactants=[Symbol(ind2reactants[ind]) for ind in 1:num_reactants]
df=DataFrame(transpose(lambda_sol)[1:end,1:num_reactants])
names!(df,reactants)
CSV.write("/data/jlbox_sensitivity_results.csv",df)
df
