#https://www.atmos-chem-phys.net/18/15743/2018/acp-18-15743-2018.pdf
#https://www.sciencedirect.com/science/article/pii/S1352231010008459
#https://www.sciencedirect.com/science/article/pii/S0187623613710622
using JlBox
using Sundials
using DataFrames
using CSV
using ArgParse

include("Generate_config.jl")

function parse_cmd()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--experiment-a", "-a"
            help = "Simulate n-th setting in experiment A"
            arg_type = Int
        "--experiment-b", "-b"
            help = "Simulate n-th setting in experiment B"
            arg_type = Int
        "--time", "-t"
            help = "Simulation period in hours"
            arg_type = Float64
            default = 1.0
        "--RH", "-r"
            help = "Specify relative humidity (out of 100 percent)"
            arg_type = Float64
        "--verbose-output"
            help = "Specify verbose output path, default to stdout"
            arg_type = String
            default = "stdout"
        "--output", "-o"
            help = "Specify path prefix for data output default to './results/'"
            arg_type = String
            default = "results/"
        "--logcpuinfo", "-i"
            help = "Log CPU info using `lscpu`"
            action = :store_true
    end
    return parse_args(s)
end


parsed_args = parse_cmd()
println(parsed_args)
iostr = parsed_args["verbose-output"]
path_pre = parsed_args["output"]
io = iostr == "stdout" ? Base.stdout : open(iostr, "w")
if parsed_args["logcpuinfo"]
    info_out = read(`lscpu`, String)
    println(io, info_out)
end
if ~xor(parsed_args["experiment-a"] isa Nothing, parsed_args["experiment-b"] isa Nothing)
    exit(-1)
elseif ~(parsed_args["experiment-a"] isa Nothing)
    time_hour = parsed_args["time"]
    name_pre = "Experiment_A"
    i = parsed_args["experiment-a"]
    config = experiment_a(time_hour, i, io, parsed_args["RH"]/100.0)
elseif ~(parsed_args["experiment-b"] isa Nothing)
    time_hour = parsed_args["time"]
    i = parsed_args["experiment-b"]
    name_pre = "Experiment_B"
    config = experiment_b(time_hour, i, io, parsed_args["RH"]/100.0)
end
println("Running Simulations for $i -th setting of $name_pre")

solverconfig = configure_aerosol_solver_sparse()
sol, reactants2ind, param_dict = JlBox.run_simulation(config, solverconfig)
if config.io != Base.stdout
    close(config.io)
end
df = JlBox.postprocess_gas(sol, reactants2ind, config)
df_SOA = JlBox.postprocess_aerosol(sol, param_dict, config)
df_size = JlBox.postprocess_aerosol_size_dist(sol, param_dict, config)
CSV.write("$(path_pre)$(name_pre)_$(i)_results.csv",df)
CSV.write("$(path_pre)$(name_pre)_$(i)_SOA.csv",df_SOA)
CSV.write("$(path_pre)$(name_pre)_$(i)_size.csv",df_size)
#julia Simulation_fullscale_manual.jl --experiment-a 1 --time 1.0
