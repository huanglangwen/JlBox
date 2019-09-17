function prepare_gas(config)
    println("Parsing Reactants")
    stoich_mtx,reactants_mtx,RO2_inds,num_eqns,num_reactants,reactants2ind=parse_reactants(config.file)
    reactants_list=mk_reactants_list(num_reactants,num_eqns,reactants_mtx)
    stoich_list=mk_reactants_list(num_reactants,num_eqns,stoich_mtx)
    @printf("num_eqns: %d, num_reactants: %d\n",num_eqns,num_reactants)

    println("Generating evaluate_rates()")
    evaluate_rates_expr=gen_evaluate_rates(config.file)
    println("Done Generation")
    rate_values=zeros(Real,num_eqns)
    #rate_prods=zeros(Float64,num_eqns)
    J=zeros(Real,62)
    #dydt=zeros(Float64,num_reactants)
    println("Performing constant folding")
    constant_folding!(evaluate_rates_expr,config.constantdict,rate_values);
    extract_constants!(evaluate_rates_expr);
    println("Evaluating evaluate_rates&loss_gain codes")
    eval(evaluate_rates_expr)
    param_dict=Dict("rate_values"=>rate_values,"J"=>J,"stoich_mtx"=>stoich_mtx,#"dydt"=>dydt,
                    "stoich_list"=>stoich_list,"reactants_list"=>reactants_list,"RO2_inds"=>RO2_inds,
                    "num_eqns"=>num_eqns,"num_reactants"=>num_reactants,"evaluate_rates!"=>evaluate_rates!,
                    "config"=>config)
    return param_dict,reactants2ind
end

function prepare_aerosol(config)
    param_dict,reactants2ind=prepare_gas(config)
    num_reactants=param_dict["num_reactants"]
    ind2reactants=Dict(reactants2ind[reac]=>reac for reac in keys(reactants2ind))
    species_names=[ind2reactants[ind] for ind=1:num_reactants]

    println("Calculating Partitioning Properties: Part1")
    pc1_dict=Pure_component1(num_reactants,species_names,config.vp_cutoff,config.temp,config.property_methods)
    
    println("Adding H2O")
    num_reactants+=1
    param_dict["num_reactants"]=num_reactants#not pc1_dict
    push!(pc1_dict["include_inds"],num_reactants)
    reactants2ind["H2O"]=num_reactants
    include_inds=pc1_dict["include_inds"]
    num_reactants_condensed=length(include_inds)
    sat_vap_water = exp(-0.58002206E4/config.temp+0.13914993E1-
        0.48640239E-1*config.temp+0.41764768E-4*(config.temp^2.0E0)-
        0.14452093E-7*(config.temp^3.0E0)+0.65459673E1*log(config.temp))#Pa
    push!(pc1_dict["y_density_array"],1000.0E0)#Append density of water to array [kg/m3]
    push!(pc1_dict["y_mw"],18.0E0)#Append mw of water to array [g/mol]
    push!(pc1_dict["Psat"],sat_vap_water*9.86923E-6)#Convert Pa to atm
    push!(pc1_dict["Delta_H"],40.66)
    Lv_water_vapour=2.5e3 # Latent heat of vapourisation of water [J/g] 
    push!(pc1_dict["Latent_heat_gas"],Lv_water_vapour)#Water vapour, taken from Paul Connolly's parcel model ACPIM

    println("Calculating Partitioning Properties: Part2")
    y_mw=pc1_dict["y_mw"]
    pc2_dict=Pure_component2(num_reactants_condensed,y_mw,config.R_gas,config.temp)
    merge!(param_dict,pc1_dict,pc2_dict)
    param_dict["num_reactants_condensed"]=num_reactants_condensed
    println("Generating initial size distribution")
    N_perbin,xs=lognormal(config.num_bins,config.total_conc,config.meansize,config.size_std,config.lowersize,config.uppersize)
    param_dict["N_perbin"]=N_perbin
    
    println("Calculating Dry Core Properties")
    y_core=(4.0/3.0)*pi*((xs*1.0e-6).^3.0) #4/3*pi*radius^3
    y_core=y_core.*config.core_density_array #mass per particle [kg]
    y_core=y_core./(config.core_mw*1.0e-3) #moles per particle, changing mw from g/mol to kg/mol
    y_core=y_core*config.NA #molecules per particle
    y_core=y_core.*N_perbin #molecules/cc representing each size range
    #Calculate a core mass based on the above information [converting from molecules/cc to micrograms/m3]    
    core_mass_array=y_core./config.NA.*config.core_mw
    println("Dry core mass = ", sum(core_mass_array)*1E12)
    param_dict["y_core"]=y_core
    param_dict["core_mass_array"]=core_mass_array

    println("Configuring initial condensed phase")
    y_cond=zeros(Float64,config.num_bins*num_reactants_condensed)
    for step=1:length(xs)
        radius=xs[step]
        water_moles=(y_core[step]*config.core_dissociation)*(config.RH/(1.0E0-config.RH))
        y_cond[step*num_reactants_condensed]=water_moles
    end
    return param_dict,reactants2ind,y_cond
end