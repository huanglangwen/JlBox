include("JlBoxModule.jl")
using .Compute:gen_simulation_gas
using TimerOutputs
using DifferentialEquations
using DifferentialEquations:CVODE_BDF

const to = TimerOutput()
prob=gen_simulation_gas()
reset_timer!(to); @timeit to "CVODE_BDF Default" begin sol=solve(prob,CVODE_BDF(),dense=false); end; show(to)
reset_timer!(to); @timeit to "CVODE_BDF GMRES" begin sol=solve(prob,CVODE_BDF(linear_solver=:GMRES),dense=false); end; show(to)
reset_timer!(to); @timeit to "Rodas4" begin sol=solve(prob,Rodas4(autodiff=false),dense=false,calck=false); end; show(to)


