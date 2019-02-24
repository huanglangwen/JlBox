include("JlBoxModule.jl")
include("Configure_aerosol.jl")
using .Compute:run_simulation_aerosol
using DataFrames
using CSV

@time sol_mtx,reactants2ind,SOA_array,num_reactants=run_simulation_aerosol()
#num_reactants=length(reactants2ind)
ind2reactants=Dict(reactants2ind[key]=>key for key in keys(reactants2ind))
reactants=[Symbol(ind2reactants[ind]) for ind in 1:num_reactants]
t_length=size(sol_mtx)[1]
t_index=linspace(0,simulation_time,t_length)
df_SOA=DataFrame(Time=t_index,SOA=SOA_array)[[:Time,:SOA]]
df=DataFrame(sol_mtx[1:end,1:num_reactants])
names!(df,reactants)
CSV.write("/data/jlbox_results.csv",df)
CSV.write("/data/jlbox_SOA.csv",df_SOA)
df_SOA
