include("JlBoxModule.jl")
using .Compute:gen_simulation_gas,to
using TimerOutputs
using DifferentialEquations
using DifferentialEquations:CVODE_BDF
using OrdinaryDiffEq
#to = TimerOutput()
prob=gen_simulation_gas()
reset_timer!(to); @timeit to "CVODE_BDF Default" begin sol=solve(prob,CVODE_BDF(),dense=false); end; show(to);println()
reset_timer!(to); @timeit to "CVODE_BDF Default2" begin sol=solve(prob,CVODE_BDF(),dense=false); end; show(to);println()
reset_timer!(to); @timeit to "CVODE_BDF GMRES" begin sol=solve(prob,CVODE_BDF(linear_solver=:GMRES),dense=false); end; show(to);println()
reset_timer!(to); @timeit to "CVODE_BDF GMRES2" begin sol=solve(prob,CVODE_BDF(linear_solver=:GMRES),dense=false); end; show(to);println()
reset_timer!(to); @timeit to "KenCarp4" begin sol=solve(prob,KenCarp4(),dense=false); end; show(to);println()