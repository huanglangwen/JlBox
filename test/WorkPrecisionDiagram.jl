include("../JlBoxModule.jl")
using .Compute:gen_simulation_gas,to
using TimerOutputs
using LinearAlgebra; BLAS.set_num_threads(8)
using OrdinaryDiffEq, Sundials, DiffEqDevTools, Plots
gr()
println("Generating Problem...")
prob,prob_jac,prob_sp=gen_simulation_gas("../data/BCR_rxn.txt","../data/BCR_pop.txt")
println("First Calibration run...")
reset_timer!(to); @timeit to "CVODE_BDF Default" begin 
    sol=solve(prob_jac,CVODE_BDF(),abstol=1/10^14,reltol=1/10^14);
end 
show(to);println()
test_sol = TestSolution(sol)

reset_timer!(to); @timeit to "Exprb32" begin sol=solve(prob_jac,Exprb32(autodiff=false)); end; show(to);println()
reset_timer!(to); @timeit to "Exprb43" begin sol=solve(prob_jac,Exprb43(autodiff=false)); end; show(to);println()

abstols = 1.0 ./ 10.0 .^ (4:7)
reltols = 1.0 ./ 10.0 .^ (1:4)

println("WorkPrecisionDiagram Run...")
setups = [Dict(:alg=>CVODE_BDF()),
          Dict(:alg=>Exprb32(autodiff=false)),
          Dict(:alg=>Exprb43(autodiff=false)),
          Dict(:alg=>Rodas4(autodiff=false)),
          Dict(:alg=>Rodas5(autodiff=false))]
names = ["CVODE_BDF" "Exprb32" "Exprb43" "Rodas4" "Rodas5"]
wp = WorkPrecisionSet(prob_jac,abstols,reltols,setups;
                      names=names,save_everystep=false,appxsol=test_sol,maxiters=Int(1e5),seconds=5)
plot(wp)