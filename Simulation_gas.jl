include("JlBoxModule.jl")
using .Compute:run_simulation_gas
using TimerOutputs
using DifferentialEquations
using DifferentialEquations:CVODE_BDF

const to = TimerOutput()
reset_timer!(to); @timeit to "Rodas4" begin sol,_=run_simulation_gas(Rodas4(autodiff=false)); end; show(to)
reset_timer!(to); @timeit to "CVODE_BDF GMRES" begin sol,_=run_simulation_gas(CVODE_BDF(linear_solver=:GMRES)); end; show(to)
reset_timer!(to); @timeit to "CVODE_BDF Default" begin sol,_=run_simulation_gas(CVODE_BDF()); end; show(to)

