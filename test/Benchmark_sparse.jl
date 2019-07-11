include("../JlBoxModule.jl")
using .Compute:gen_simulation_gas,to
using TimerOutputs
#using Sundials:CVODE_BDF
using OrdinaryDiffEq

prob,prob_jac,prob_sp=gen_simulation_gas("../data/BCR_rxn.txt","../data/BCR_pop.txt")
reset_timer!(to); @timeit to "KenCarp4 FD" begin sol=solve(prob,KenCarp4(autodiff=false),dense=false); end; show(to);println()
reset_timer!(to); @timeit to "KenCarp4 FD sparse" begin sol=solve(prob_sp,KenCarp4(autodiff=false),dense=false); end; show(to);println()
reset_timer!(to); @timeit to "KenCarp4 AD" begin sol=solve(prob,KenCarp4(),dense=false); end; show(to);println()
reset_timer!(to); @timeit to "KenCarp4 AD sparse" begin sol=solve(prob_sp,KenCarp4(),dense=false); end; show(to);println()
reset_timer!(to); @timeit to "KenCarp4 jac" begin sol=solve(prob_jac,KenCarp4(autodiff=false),dense=false); end; show(to);println()