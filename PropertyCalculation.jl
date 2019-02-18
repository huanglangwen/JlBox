#Direct Interpreted from https://github.com/loftytopping/PyBox/blob/master/Aerosol/Property_calculation.py
using PyCall
using LightXML
unshift!(PyVector(pyimport("sys")["path"]),"../UManSysProp_public")
@pyimport umansysprop.boiling_points as boiling_points
@pyimport umansysprop.vapour_pressures as vapour_pressures
@pyimport umansysprop.critical_properties as critical_properties
@pyimport umansysprop.liquid_densities as liquid_densities
@pyimport umansysprop.partition_models as partition_models
@pyimport umansysprop.groups as groups
@pyimport umansysprop.activity_coefficient_models_dev as aiomfac
@pyimport umansysprop.forms as forms #need forms.CoreAbundanceField (class)
CoreAbundanceField=forms.CoreAbundanceField
@pyimport pybel

function readSMILESdict()
    species2SMILESdict=Dict{String,String}()
    mcmdoc=parse_file("MCM.xml")
    spnode=find_element(root(mcmdoc),"species_defs")
    for species in child_nodes(spnode)
        if is_elementnode(species)
            sp_element=XMLElement(species)
            name=attribute(sp_element,"species_name")
            if has_children(sp_element)
                SMILES=strip(content(sp_element))
                species2SMILESdict[name]=SMILES
            end
        end
    end
    return species2SMILESdict
end

function SMILES2Pybel(smi_str)
    return pybel.readstring("smi",smi_str)
end

function compoundProperty(pybelobj::PyObject,temperature::Integer,methodfuncs::Dict)
    boiling_point,vapour_pressure,critical_property,liquid_density=[methodfuncs[i] for i in ["bp","vp","critical","density"]]
    b1=boiling_points.nannoolal(pybelobj)
    density=liquid_density(pybelobj, temperature, pycall(critical_property,PyObject,pybelobj, b1))*1.0E3
    mw=pybelobj[:molwt]
    groups_dict=groups.composition(pybelobj)
    o_c=groups_dict["O"]/groups_dict["C"]
    h_c=groups_dict["H"]/groups_dict["C"]
    b=boiling_point(pybelobj)
    sat_vp=vapour_pressure(pybelobj, temperature,b)
    return density,mw,o_c,h_c,sat_vp
end

function Pure_component1(num_species::Integer,species_names::Array{String,1},vp_cutoff::Real,temperature::Number,methods::Dict)
    species2inds=Dict(species_names[i]=>i for i in 1:num_species)
    bp_method,vp_method,critical_method,density_method=[methods[i] for i in ["bp","vp","critical","density"]]
    boiling_point = Dict(
        "joback_and_reid"=> boiling_points.joback_and_reid,
        "stein_and_brown"=> boiling_points.stein_and_brown,
        "nannoolal"=>       boiling_points.nannoolal,
    )[bp_method]
    vapour_pressure = Dict(
        "nannoolal"=>            vapour_pressures.nannoolal,
        "myrdal_and_yalkowsky"=> vapour_pressures.myrdal_and_yalkowsky,
        # Evaporation doesn"t use boiling point
        "evaporation"=> (c, t, b)-> vapour_pressures.evaporation(c, t),
    )[vp_method]
    critical_property = Dict(
        "nannoolal"=>          critical_properties.nannoolal,
        "joback_and_reid"=>    critical_properties.joback_and_reid,
    )[critical_method]
    liquid_density = Dict(
        "girolami"=>      (c, t, p) -> liquid_densities.girolami(c),
        "schroeder"=>     liquid_densities.schroeder,
        "le_bas"=>        liquid_densities.le_bas,
        "tyn_and_calus"=> liquid_densities.tyn_and_calus,
    )[density_method]
    methodfuncs=Dict(
        "bp"=>boiling_point,
        "vp"=>vapour_pressure,
        "critical"=>critical_property,
        "density"=>liquid_density
    )
    y_density_array=zeros(Float64,num_species)+1000
    y_mw=zeros(Float64,num_species)+200
    o_cs=zeros(Float64,num_species)
    h_cs=zeros(Float64,num_species)
    sat_vps=zeros(Float64,num_species)+100
    #ignore_index_mask=zeros(Bool,num_species)#All falses
    include_inds::Array{Integer,1}=[]
    Delta_H=zeros(Float64,num_species)#Ignored
    Latent_heat_gas=zeros(Float64,num_species)#Ignored

    println("Reading SMILES definitions")
    species2SMILESdict=readSMILESdict()

    println("Calculating component properties using UManSysProp")
    for species_ind in 1:num_species
        species_name=species_names[species_ind]
        if haskey(species2SMILESdict,species_name)
            pybelobj=SMILES2Pybel(species2SMILESdict[species_name])
            density,mw,o_c,h_c,sat_vp=compoundProperty(pybelobj,temperature,methodfuncs)
            if (typeof(density)<:Real)&(sat_vp<=vp_cutoff)
                push!(include_inds,species_ind)
                y_density_array[species_ind]=density
                y_mw[species_ind]=mw
                o_cs[species_ind]=o_c
                h_cs[species_ind]=h_c
                sat_vps[species_ind]=sat_vp
            end
        end
    end
    return_dict=Dict(
        "y_density_array"=>y_density_array[include_inds],
        "y_mw"=>y_mw[include_inds],
        "o_c"=>o_cs[include_inds],
        "h_c"=>h_cs[include_inds],
        #"sat_vp"=>sat_vps[include_inds],
        "Psat"=>exp10.(sat_vps[include_inds]),
        "Delta_H"=>Delta_H[include_inds],
        "Latent_heat_gas"=>Latent_heat_gas[include_inds],
        "include_inds"=>include_inds
        #"ignore_index_mask"=>ignore_index_mask
    )
end

function Pure_component2(num_species_condensed::Integer,y_mw::Array{Float64,1},R_gas::Real,temperature::Real)::Dict{String,Array{Real,1}}
    alpha_d_org=zeros(Float64,num_species_condensed)+0.1
    DStar_org=1.9E0*(y_mw.^(-2.0E0/3.0E0))
    mean_them_vel=sqrt.((8.0E0*R_gas*temp)./(pi*y_mw*1.0E-3))
    gamma_gas=((3.0E0*DStar_org)./(mean_them_vel*1.0E2))*1.0E-2

    return_dict=Dict(
        "alpha_d_org"=>alpha_d_org,
        "DStar_org"=>DStar_org,
        "mean_them_vel"=>mean_them_vel,
        "gamma_gas"=>gamma_gas
    )
    return return_dict
end
