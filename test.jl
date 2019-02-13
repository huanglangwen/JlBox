include("PropertyCalculation.jl")
function test_Pure_component1()
    sp2SMILES=readSMILESdict()
    species_names=collect(keys(sp2SMILES))
    num_species=length(species_names)
    temperature=298
    methods=Dict("bp"=>"nannoolal","vp"=>"nannoolal","critical"=>"nannoolal","density"=>"schroeder")
    Pure_component1(num_species,species_names,temperature,methods)
end
