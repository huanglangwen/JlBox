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
