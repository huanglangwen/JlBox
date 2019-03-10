include("PropertyCalculation.jl")
function test_Pure_component1()
    sp2SMILES=readSMILESdict()
    species_names=collect(keys(sp2SMILES))
    num_species=length(species_names)
    temperature=298
    vp_cutoff=-6.0
    methods=Dict("bp"=>"nannoolal","vp"=>"nannoolal","critical"=>"nannoolal","density"=>"schroeder")
    Pure_component1(num_species,species_names,vp_cutoff,temperature,methods)
end

include("JlBoxModule.jl")
using .Compute:read_configure!,prepare_aerosol
function test_aerosol_initial()
    include("Configure_aerosol.jl")
    read_configure!("Configure_aerosol.jl")
    param_dict,reactants2ind,y_cond=prepare_aerosol()
    num_reactants,num_reactants_condensed=[param_dict[i] for i in ["num_reactants","num_reactants_condensed"]]
    dy_dt_gas_matrix=zeros(Float64,(num_reactants,num_bins))
    dy_dt=zeros(Float64,num_reactants+num_reactants_condensed*num_bins)
    param_dict["dy_dt_gas_matrix"]=dy_dt_gas_matrix
    param_dict["dydt"]=dy_dt
    param_dict["Current_iter"]=0
    y_init=zeros(Float64,num_reactants+num_reactants_condensed*num_bins)
    for (k,v) in reactants_initial_dict
        y_init[reactants2ind[k]]=v*Cfactor#pbb to molcules/cc
    end
    return param_dict,reactants2ind
end

function test_include_inds()
    param_dict,reactants2ind=test_aerosol_initial()
    ind2reactants=Dict(reactants2ind[key]=>key for key in keys(reactants2ind))
    num_reactants=param_dict["num_reactants"]
    reactants=[ind2reactants[ind] for ind in 1:num_reactants]
    include_inds=param_dict["include_inds"]

    sp2SMILES=readSMILESdict()
    for i in [47,48,49]
        reac_name=reactants[i]
        println("reactant_name:",reac_name)
        if haskey(sp2SMILES,reac_name)
            reac_smi=sp2SMILES[reac_name]
            reac_pyobj=SMILES2Pybel(reac_smi)
            println("SMILES:",reac_smi)
            
            property_methods=Dict("bp"=>boiling_points.joback_and_reid,
                                  "vp"=>vapour_pressures.nannoolal,
                                  "critical"=>critical_properties.nannoolal,
                                  "density"=>(c, t, p) -> liquid_densities.girolami(c))
            density,mw,o_c,h_c,sat_vp=compoundProperty(reac_pyobj,temp,property_methods)
            println(density," ,",mw," ,",o_c," ,",h_c," ,",sat_vp)
        else
            println("not in MCM.xml")
        end
    end
end

using DataFrames
using CSV
function test_properties()
    param_dict,reactants2ind=test_aerosol_initial()
    props_keys=["y_density_array","y_mw","Psat",
                "Delta_H","Latent_heat_gas","include_inds",
                "alpha_d_org","DStar_org","gamma_gas"]
    props_symbols=[Symbol(k) for k in props_keys]
    props_dict=filter((k,v)->k in props_keys,param_dict)
    props_dict["sat_vp"]=log10.(props_dict["Psat"])
    props_df=DataFrame(props_dict)[props_symbols]
    CSV.write("/data/jlbox_props.csv",props_df)
    nothing
end

using .Compute:read_configure!,prepare_gas,loss_gain!
using .Jacobian:loss_gain_jac!
function test_jacobian()
    include("Configure_gas.jl")
    read_configure!("Configure_gas.jl")
    param_dict,reactants2ind,evaluate_rates_expr=prepare_gas()
    eval(evaluate_rates_expr)
    rate_values,J,stoich_mtx,stoich_list,reactants_list,RO2_inds,num_eqns,num_reactants=
        [param_dict[ind] for ind in 
            ["rate_values","J","stoich_mtx","stoich_list","reactants_list","RO2_inds",
             "num_eqns","num_reactants"]
        ]
    ind2reactants=Dict(reactants2ind[key]=>key for key in keys(reactants2ind))
    reactants_names=[ind2reactants[ind] for ind in 1:num_reactants]
    reactants_initial=zeros(Float64,num_reactants)
    dydt=zeros(Float64,num_reactants)
    for (k,v) in reactants_initial_dict
        reactants_initial[reactants2ind[k]]=v*Cfactor#pbb to molcules/cc
    end
    time_of_day_seconds=start_time+0
    RO2=sum(reactants_initial[RO2_inds])
    evaluate_rates!(time_of_day_seconds,RO2,H2O,temp,rate_values,J)# =>ratevalues
    dydt_raw=deepcopy(loss_gain!(num_reactants,num_eqns,reactants_initial,stoich_mtx,stoich_list,reactants_list,rate_values,dydt))
    delta=1E-10
    invdelta=1E10
    lossgain_jac_mtx=zeros(num_reactants,num_reactants)#num_output(dydt)*num_input(y)
    lossgain_jac_mtx2=zeros(num_reactants,num_reactants)
    inc_array=zeros(Float64,num_reactants)
    for reactant_ind in 1:num_reactants
        if reactant_ind>=2
            inc_array[reactant_ind-1]=0
        end
        inc_array[reactant_ind]=delta
        loss_gain!(num_reactants,num_eqns,reactants_initial.+inc_array,stoich_mtx,stoich_list,reactants_list,rate_values,dydt)
        lossgain_jac_mtx[:,reactant_ind]=(dydt.-dydt_raw).*invdelta
    end
    loss_gain_jac!(num_reactants,num_eqns,reactants_initial,stoich_mtx,stoich_list,reactants_list,rate_values,lossgain_jac_mtx2)
    df1=DataFrame(lossgain_jac_mtx)
    df2=DataFrame(lossgain_jac_mtx2)
    names!(df1,[Symbol(reac) for reac in reactants_names])
    names!(df2,[Symbol(reac) for reac in reactants_names])
    CSV.write("/data/lossgain_jac1.csv",df1)
    CSV.write("/data/lossgain_jac2.csv",df2)
    df1,df2
end

function test_aerosol_jacobian()
    include("Configure_aerosol.jl")
    read_configure!("Configure_aerosol.jl")
    param_dict,reactants2ind,y_cond,evaluate_rates_expr=prepare_aerosol()
    nothing
end