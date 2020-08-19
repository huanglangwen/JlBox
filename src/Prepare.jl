function get_sparsity_gas(param_dict,reactants2ind)
    rate_values,J,stoich_mtx,stoich_list,reactants_list,RO2_inds,num_eqns,num_reactants=
        [param_dict[ind] for ind in 
            ["rate_values","J","stoich_mtx","stoich_list","reactants_list","RO2_inds",
             "num_eqns","num_reactants"]
        ]
    config=param_dict["config"]
    evaluate_rates_fun=param_dict["evaluate_rates!"]
    time_of_day_seconds=config.start_time
    RO2=1e5
    Base.invokelatest(evaluate_rates_fun,time_of_day_seconds,RO2,config.H2O,config.temp,rate_values,J)# =>ratevalues
    jac_prototype=zeros(num_reactants,num_reactants)
    y_init=ones(Float64,num_reactants)*1.5
    for (k,v) in config.reactants_initial_dict
        y_init[reactants2ind[k]]=v*config.Cfactor#pbb to molcules/cc
    end
    loss_gain_jac!(num_reactants,num_eqns,y_init,stoich_mtx,stoich_list,reactants_list,rate_values,jac_prototype)
    sparse(jac_prototype)
end

function get_sparsity_aerosol(solverconfig::SolverConfig,param_dict,reactants2ind,y_cond)
    config,num_reactants,num_reactants_condensed=[param_dict[i] for i in ["config","num_reactants","num_reactants_condensed"]]
    len_y=num_reactants+num_reactants_condensed*config.num_bins
    y_init=ones(Float64,len_y)*1.5
    for (k,v) in config.reactants_initial_dict
        y_init[reactants2ind[k]]=v*config.Cfactor#pbb to molcules/cc
    end
    y_init[num_reactants+1:num_reactants+config.num_bins*num_reactants_condensed]=y_cond[1:config.num_bins*num_reactants_condensed]
    jac! = select_jacobian(solverconfig.diff_method,len_y)
    jac_prototype=zeros(len_y,len_y)
    jac!(jac_prototype,y_init,param_dict,0.0)
    sparse(jac_prototype)
end

function prepare_gas(config::JlBoxConfig)
    io = config.io
    @debug "Parsing Reactants"
    stoich_mtx,reactants_mtx,RO2_inds,num_eqns,num_reactants,reactants2ind=parse_reactants(config.file)
    reactants_list=mk_reactants_list(num_reactants,num_eqns,reactants_mtx)
    stoich_list=mk_reactants_list(num_reactants,num_eqns,stoich_mtx)
    println(io, "num_eqns: $(num_eqns), num_reactants: $(num_reactants)")

    @debug "Generating evaluate_rates()"
    evaluate_rates_expr=gen_evaluate_rates(config)
    @debug "Done Generation"
    rate_values=zeros(Real,num_eqns)
    J=zeros(Real,62)
    @debug "Performing constant folding"
    constant_folding!(evaluate_rates_expr,config.constant_dict,rate_values);
    extract_constants!(evaluate_rates_expr);
    @debug "Evaluating evaluate_rates&loss_gain codes"
    eval(evaluate_rates_expr)
    param_dict=Dict("rate_values"=>rate_values,"J"=>J,"stoich_mtx"=>stoich_mtx,#"dydt"=>dydt,
                    "stoich_list"=>stoich_list,"reactants_list"=>reactants_list,"RO2_inds"=>RO2_inds,
                    "num_eqns"=>num_eqns,"num_reactants"=>num_reactants,"evaluate_rates!"=>evaluate_rates!,
                    "config"=>config)
    return param_dict,reactants2ind
end

function prepare_aerosol(config::JlBoxConfig,param_dict::Dict,reactants2ind)
    io = config.io
    num_reactants=param_dict["num_reactants"]
    ind2reactants=Dict(reactants2ind[reac]=>reac for reac in keys(reactants2ind))
    species_names=[ind2reactants[ind] for ind=1:num_reactants]

    @debug "Calculating Partitioning Properties: Part1"
    pc1_dict=Pure_component1(num_reactants,species_names,config.vp_cutoff,config.temp,config.property_methods)
    
    @debug "Adding H2O"
    num_reactants+=1
    param_dict["num_reactants"]=num_reactants#not pc1_dict
    push!(pc1_dict["include_inds"],num_reactants)
    reactants2ind["H2O"]=num_reactants
    include_inds=pc1_dict["include_inds"]
    num_reactants_condensed=length(include_inds)
    println(io, "num_reactants_condensed: $(num_reactants_condensed)")
    sat_vap_water = exp(-0.58002206E4/config.temp+0.13914993E1-
        0.48640239E-1*config.temp+0.41764768E-4*(config.temp^2.0E0)-
        0.14452093E-7*(config.temp^3.0E0)+0.65459673E1*log(config.temp))#Pa
    push!(pc1_dict["y_density_array"],1000.0E0)#Append density of water to array [kg/m3]
    push!(pc1_dict["y_mw"],18.0E0)#Append mw of water to array [g/mol]
    push!(pc1_dict["Psat"],sat_vap_water*9.86923E-6)#Convert Pa to atm
    push!(pc1_dict["Delta_H"],40.66)
    Lv_water_vapour=2.5e3 # Latent heat of vapourisation of water [J/g] 
    push!(pc1_dict["Latent_heat_gas"],Lv_water_vapour)#Water vapour, taken from Paul Connolly's parcel model ACPIM

    @debug "Calculating Partitioning Properties: Part2"
    y_mw=pc1_dict["y_mw"]
    pc2_dict=Pure_component2(num_reactants_condensed,y_mw,config.temp)
    merge!(param_dict,pc1_dict,pc2_dict)
    param_dict["num_reactants_condensed"]=num_reactants_condensed
    @debug "Generating initial size distribution"
    N_perbin,xs=lognormal(config.num_bins,config.total_conc,config.meansize,config.size_std,config.lowersize,config.uppersize)
    param_dict["N_perbin"]=N_perbin
    
    @debug "Calculating Dry Core Properties"
    y_core=(4.0/3.0)*pi*((xs*1.0e-6).^3.0) #4/3*pi*radius^3
    y_core=y_core.*config.core_density_array #mass per particle [kg]
    y_core=y_core./(config.core_mw*1.0e-3) #moles per particle, changing mw from g/mol to kg/mol
    y_core=y_core*NA #molecules per particle
    y_core=y_core.*N_perbin #molecules/cc representing each size range
    #Calculate a core mass based on the above information [converting from molecules/cc to micrograms/m3]    
    core_mass_array=y_core./NA.*config.core_mw
    println(io, "Dry core mass = ", sum(core_mass_array)*1E12)
    param_dict["y_core"]=y_core
    param_dict["core_mass_array"]=core_mass_array

    @debug "Configuring initial condensed phase"
    y_cond=zeros(Float64,config.num_bins*num_reactants_condensed)
    for step=1:length(xs)
        radius=xs[step]
        water_moles=(y_core[step]*config.core_dissociation)*(config.RH/(1.0E0-config.RH))
        y_cond[step*num_reactants_condensed]=water_moles
    end
    dy_dt_gas_matrix=zeros(Real,(num_reactants,config.num_bins))
    param_dict["dy_dt_gas_matrix"]=dy_dt_gas_matrix
    #set initial value
    len_y=num_reactants+num_reactants_condensed*config.num_bins
    y_init=zeros(Float64,len_y)
    for (k,v) in config.reactants_initial_dict
        y_init[reactants2ind[k]]=v*config.Cfactor#pbb to molcules/cc
    end
    y_init[num_reactants+1:num_reactants+config.num_bins*num_reactants_condensed]=y_cond[1:config.num_bins*num_reactants_condensed]
    return param_dict,reactants2ind,y_init
end

function prepare(config::JlBoxConfig, solverconfig::SolverConfig)
    if config isa GasConfig
        param_dict, reactants2ind = prepare_gas(config)
        param_dict["ShowIterPeriod"] = 100
        #Set initial value
        y_init = zeros(Float64,param_dict["num_reactants"])
        for (k,v) in config.reactants_initial_dict
            y_init[reactants2ind[k]]=v*config.Cfactor#pbb to molcules/cc
        end
    elseif config isa AerosolConfig
        param_dict, reactants2ind = prepare_gas(config)
        param_dict, reactants2ind, y_init = prepare_aerosol(config, param_dict, reactants2ind)
        param_dict["ShowIterPeriod"] = 10
    end
    param_dict["Current_iter"] = 0
    len_y = length(y_init)
    jac_prototype = solverconfig.sparse ? spzeros(len_y, len_y) : nothing
    param_dict["sparsity"] = jac_prototype
    return param_dict, reactants2ind, y_init
end
