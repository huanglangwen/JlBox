using JlBox
using Sundials
import MPI
using DataFrames
using CSV

MPI.Init()

comm = MPI.COMM_WORLD
MPI.Barrier(comm)

include("Generate_config.jl")

rank = MPI.Comm_rank(comm)
nrank = MPI.Comm_size(comm)

time_hour = 1.0
#configs = vcat(experiment_a(time_hour), experiment_b(time_hour))
if nrank == 7
    configs = experiment_b(time_hour)
    name = "Experiment_B_rank" * rank
elseif nrank == 9
    configs = experiment_a(time_hour)
    name = "Experiment_A_rank" * rank
else
    exit(-1)
end

config = configs[rank]
solverconfig = configure_aerosol_solver_sparse()
#rets = [remotecall((i-1)%nranks, JlBox.run_simulation, configs[i], solverconfig) for i in 1:16]
sol, reactants2ind, param_dict = JlBox.run_simulation(config, solverconfig)
df = JlBox.postprocess_gas(sol, reactants2ind, config)
df_SOA = JlBox.postprocess_aerosol(sol, param_dict, config)
CSV.write("result/"*name*"_results.csv",df)
CSV.write("result/"*name*"_SOA.csv",df_SOA)

# using MPI; MPI.install_mpiexecjl()
# mpiexecjl --output-filename result/log -n 7 julia Simulation_fullscale.jl