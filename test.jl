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
include("Configure_aerosol.jl")
using Compute:read_configure!,prepare_aerosol
function test_aerosol_initial()
    read_configure!("Configure_aerosol.jl")
    param_dict,reactants2ind=prepare_aerosol()
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

function test_properties()
    param_dict,reactants2ind=test_aerosol_initial()
    
end