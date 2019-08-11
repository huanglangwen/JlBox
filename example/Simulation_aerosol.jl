using JlBox
using DataFrames
#using CSV

@time sol,reactants2ind,SOA_array,num_reactants,_=run_simulation_aerosol(use_jacobian=false)
#num_reactants=length(reactants2ind)
sol_mtx=transpose(sol)
ind2reactants=Dict(reactants2ind[key]=>key for key in keys(reactants2ind))
reactants=[Symbol(ind2reactants[ind]) for ind in 1:num_reactants]
t_length=size(sol_mtx)[1]
t_index=range(0,stop=simulation_time,length=t_length)
df_SOA=DataFrame(Time=t_index,SOA=SOA_array)[[:Time,:SOA]]
df=DataFrame(sol_mtx[1:end,1:num_reactants])
names!(df,reactants)
#CSV.write("/data/jlbox_results.csv",df)
#CSV.write("/data/jlbox_SOA.csv",df_SOA)
df_SOA
