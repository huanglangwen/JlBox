using DataFrames
function postprocess_gas(sol, reactants2ind, config::JlBoxConfig)
    num_reactants = length(reactants2ind)
    ind2reactants = Dict(reactants2ind[key]=>key for key in keys(reactants2ind))
    reactants = [ind2reactants[ind] for ind in 1:num_reactants]
    df = DataFrames.DataFrame(transpose(sol)[1:end,1:num_reactants], :auto)
    t_length = size(df, 1)
    t_index = range(0, stop = config.simulation_time, length = t_length)
    DataFrames.rename!(df,[Symbol(reac) for reac in reactants])
    DataFrames.insertcols!(df, 1, :Time => t_index)
    return df
end

function postprocess_aerosol(sol, param_dict, config::JlBoxConfig)
    sol_mtx = transpose(sol)
    num_reactants = param_dict["num_reactants"]
    num_reactants_condensed = param_dict["num_reactants_condensed"]
    num_bins = convert(Int, (size(sol_mtx)[2] - num_reactants)/num_reactants_condensed)
    aerosol_mtx = sol_mtx[1:end,num_reactants+1:num_reactants+num_bins*num_reactants_condensed]
    t_length = size(aerosol_mtx)[1]
    t_index = range(0, stop = config.simulation_time, length = t_length)
    mw_array = param_dict["y_mw"]
    core_mass_array = param_dict["core_mass_array"]
    SOA_array = [sum((sum(reshape(aerosol_mtx[i,1:end],(num_reactants_condensed,num_bins))
                               ,dims=2).*mw_array./NA)[1:end-1]#exclude H2O at the end
                  ) for i in 1:t_length]*1E12
    mass_mtx = zeros((t_length, num_bins))
    for i in 1:t_length
        mass_mtx[i,:] = transpose(mw_array)*reshape(aerosol_mtx[i,1:end],(num_reactants_condensed,num_bins))./NA
        mass_mtx[i,:] += core_mass_array
        mass_mtx[i,:] .*= 1E12 #ug/m3
    end
    df_SOA = DataFrames.DataFrame(mass_mtx, :auto)
    DataFrames.rename!(df_SOA,["Bin_$(i)" for i in 1:num_bins])
    DataFrames.insertcols!(df_SOA, 1, :Time => t_index, :SOA => SOA_array)
    return df_SOA
end

function compute_size(num_bins, num_reactants, num_reactants_condensed, y, mw_array, N_perbin, density_input, core_density_array, core_mass_array)
    mass_array=zeros(eltype(y),num_reactants_condensed+1)
    density_array=zeros(eltype(density_input),num_reactants_condensed+1)
    size_array=zeros(num_bins)
    for size_step=1:num_bins
        start_ind=num_reactants+1+((size_step-1)*num_reactants_condensed)
        stop_ind=num_reactants+(size_step*num_reactants_condensed)
        temp_array=y[start_ind:stop_ind]
        mass_array[1:num_reactants_condensed]=temp_array.*mw_array./NA
        mass_array[num_reactants_condensed+1]=core_mass_array[size_step]
        density_array[1:num_reactants_condensed]=density_input[1:num_reactants_condensed]
        density_array[num_reactants_condensed+1]=core_density_array[size_step]
        
        #aw_array[size_step]=temp_array[num_reactants_condensed]/total_moles
        total_mass=sum(mass_array)
        total_mass= total_mass > 0 ? total_mass : core_mass_array[size_step]#fix negative size_array^3
        mass_fractions_array=mass_array./total_mass

        #density=1.0/(sum(mass_fractions_array./density_array))
        invdensity=sum(mass_fractions_array./density_array)
        invdensity= invdensity > 0 ? invdensity : 1/core_density_array[size_step]
        
        size_array[size_step]=((3.0*((total_mass*1.0E3)/(N_perbin[size_step]*1.0E6)))*invdensity/(4.0*pi))^(1.0/3.0)
    end
    return size_array
end

function postprocess_aerosol_size_dist(sol, param_dict, config::JlBoxConfig)
    sol_mtx = transpose(sol)
    num_reactants = param_dict["num_reactants"]
    num_reactants_condensed = param_dict["num_reactants_condensed"]
    num_bins = convert(Int, (size(sol_mtx)[2] - num_reactants)/num_reactants_condensed)
    aerosol_mtx = sol_mtx[1:end,num_reactants+1:num_reactants+num_bins*num_reactants_condensed]
    t_length = size(aerosol_mtx)[1]
    t_index = range(0, stop = config.simulation_time, length = t_length)
    mw_array = param_dict["y_mw"]
    N_perbin = param_dict["N_perbin"]
    density_input = param_dict["y_density_array"]
    core_mass_array = param_dict["core_mass_array"]
    size_mtx = zeros((t_length, num_bins))
    for i in 1:t_length
        y = sol_mtx[i,:]
        size_mtx[i,:] = compute_size(num_bins, num_reactants, num_reactants_condensed, y, mw_array, N_perbin, density_input, config.core_density_array, core_mass_array)
    end
    df_size = DataFrames.DataFrame(size_mtx)
    DataFrames.rename!(df_size,["Bin_$(i)" for i in 1:num_bins])
    DataFrames.insertcols!(df_size, 1, :Time => t_index)
    return df_size
end