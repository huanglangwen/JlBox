module Optimize
export mk_stoich_list
include("Optimize.jl")
end

module Parse_eqn
export parse_reactants
include("Parse_eqn.jl")
end

module Jacobian
export gas_jac!
include("Jacobian.jl")
end

module Compute
export run_simulation_gas,run_simulation_aerosol
include("Compute.jl")
end