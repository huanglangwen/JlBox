include("JlBoxModule.jl")
using .Compute:gen_simulation_gas
using TimerOutputs
using DifferentialEquations
using DifferentialEquations:CVODE_BDF

const to = TimerOutput()
prob=gen_simulation_gas()
reset_timer!(to); @timeit to "ESDIRK" begin sol=solve(prob,ESDIRK54I8L2SA(autodiff=false, diff_type=Val{:forward}, new_jac_conv_bound=Inf), dt=0.16); end; show(to)
reset_timer!(to); @timeit to "CVODE_BDF Default" begin sol=solve(prob,CVODE_BDF(),dense=false); end; show(to)
reset_timer!(to); @timeit to "CVODE_BDF GMRES" begin sol=solve(prob,CVODE_BDF(linear_solver=:GMRES),dense=false); end; show(to)