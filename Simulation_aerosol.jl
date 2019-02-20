include("JlBoxModule.jl")
include("Configure_aerosol.jl")
using Compute:run_simulation_aerosol
using DataFrames
using CSV

@time sol,reactants2ind,num_reactants,num_reactants_condensed=run_simulation_aerosol()
#num_reactants=length(reactants2ind)
ind2reactants=Dict(reactants2ind[key]=>key for key in keys(reactants2ind))
reactants=[Symbol(ind2reactants[ind]) for ind in 1:num_reactants]
sol_mtx=transpose(sol)
aerosol_mtx=sol_mtx[1:end,num_reactants+1:num_reactants+num_bins*num_reactants_condensed]
t_length=size(aerosol_mtx)[1]
SOA_array=[sum(sum(reshape(aerosol_mtx[i,1:end],(num_reactants_condensed,num_bins)),2)[1:end-1]) for i in 1:t_length]
t_index=linspace(0,3600,t_length)
df_SOA=DataFrame(Time=t_index,SOA_array)
df=DataFrame(sol_mtx[1:end,1:num_reactants])
names!(df,reactants)
CSV.write("/data/jlbox_results.csv",df)
df
