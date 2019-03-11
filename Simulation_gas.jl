include("JlBoxModule.jl")
using .Compute:gen_simulation_gas,to
using TimerOutputs
using DifferentialEquations
using DifferentialEquations:CVODE_BDF
using OrdinaryDiffEq
#to = TimerOutput()
prob,prob_jac=gen_simulation_gas()
reset_timer!(to); @timeit to "CVODE_BDF Default" begin sol=solve(prob,CVODE_BDF(),dense=false,abstol=1e-8, reltol=1e-8); end; show(to);println()
reset_timer!(to); @timeit to "CVODE_BDF Default with Jac" begin sol=solve(prob_jac,CVODE_BDF(),dense=false,abstol=1e-8, reltol=1e-8); end; show(to);println()
reset_timer!(to); @timeit to "CVODE_BDF GMRES" begin sol=solve(prob,CVODE_BDF(linear_solver=:GMRES),dense=false,abstol=1e-8, reltol=1e-8); end; show(to);println()
reset_timer!(to); @timeit to "ROCK2" begin sol=solve(prob,ROCK2(),dense=false); end; show(to);println()
reset_timer!(to); @timeit to "ROCK4" begin sol=solve(prob,ROCK4(),dense=false); end; show(to);println()