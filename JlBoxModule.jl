module Optimize
export constant_folding!,extract_constants!,generate_loss_gain,mk_stoich_list
include("Optimize.jl")
end

module Parse_eqn
export parse_reactants,gen_evaluate_rates
include("Parse_eqn.jl")
end

module Compute
export run_simulation_gas,run_simulation_aerosol
include("Compute.jl")
end