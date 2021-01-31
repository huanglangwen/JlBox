#https://www.atmos-chem-phys.net/18/15743/2018/acp-18-15743-2018.pdf
using JlBox
using Sundials
using Distributed
using DataFrames
using CSV

include("Generate_config.jl")

@everywhere function execute_rank_n(config, solverconfig, rank, name_pre)
    sol, reactants2ind, param_dict = JlBox.run_simulation(config, solverconfig)
    if config.io != Base.stdout
        close(config.io)
    end
    df = JlBox.postprocess_gas(sol, reactants2ind, config)
    df_SOA = JlBox.postprocess_aerosol(sol, param_dict, config)
    CSV.write("results/$(name_pre)$(rank)_results.csv",df)
    CSV.write("results/$(name_pre)$(rank)_SOA.csv",df_SOA)
end


nrank = nworkers()
println("Running Simulations with $(rank) workers")
time_hour = 1.0
configs = vcat(experiment_a(time_hour), experiment_b(time_hour))

names = [i<=9 ? "Experiment_A_rank" : "Experiment_B_rank" for i in 1:16]

solverconfig = configure_aerosol_solver_sparse()
rets = [remotecall(execute_rank_n, (i-1)%nrank + 1, configs[i], solverconfig, i, names[i]) for i in 1:16]
for ret in rets
    fetch(ret)
end
#julia -p 7 Simulation_fullscale_distributed.jl