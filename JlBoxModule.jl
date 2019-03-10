module Optimize
export mk_stoich_list
include("Optimize.jl")
end

module Parse_eqn
export parse_reactants
include("Parse_eqn.jl")
end

module Jacobian
export loss_gain_jac!
include("Jacobian.jl")
end

module Compute
include("Compute.jl")
end