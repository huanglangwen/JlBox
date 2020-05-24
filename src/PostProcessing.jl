using DataFrames
function postprocess_gas(sol, reactants2ind)
    num_reactants = length(reactants2ind)
    ind2reactants = Dict(reactants2ind[key]=>key for key in keys(reactants2ind))
    reactants = [ind2reactants[ind] for ind in 1:num_reactants]
    df = DataFrames.DataFrame(transpose(sol)[1:end,1:num_reactants])
    DataFrames.rename!(df,[Symbol(reac) for reac in reactants])
    return df
end

function postprocess_aerosol(sol, param_dict, simulation_time)
    sol_mtx = transpose(sol)
    num_reactants = param_dict["num_reactants"]
    num_reactants_condensed = param_dict["num_reactants_condensed"]
    num_bins = convert(Int, (size(sol_mtx)[2] - num_reactants)/num_reactants_condensed)
    aerosol_mtx = sol_mtx[1:end,num_reactants+1:num_reactants+num_bins*num_reactants_condensed]
    t_length = size(aerosol_mtx)[1]
    t_index = range(0, stop = simulation_time, length = t_length)
    mw_array = param_dict["y_mw"]
    SOA_array = [sum((sum(reshape(aerosol_mtx[i,1:end],(num_reactants_condensed,num_bins))
                               ,dims=2).*mw_array./NA)[1:end-1]#exclude H2O at the end
                  ) for i in 1:t_length]*1E12
    df_SOA = DataFrames.DataFrame(Time=t_index,SOA=SOA_array)[:,[:Time,:SOA]]
    return df_SOA
end