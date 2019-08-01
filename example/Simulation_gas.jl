include("JlBoxModule.jl")
using .Compute:run_simulation_gas
using DataFrames
using CSV

#Profile.init(n = 10^7, delay = 5.)
sol,reactants2ind=run_simulation_gas(use_jacobian=true)
num_reactants=length(reactants2ind)
ind2reactants=Dict(reactants2ind[key]=>key for key in keys(reactants2ind))
reactants=[ind2reactants[ind] for ind in 1:num_reactants]
df=DataFrame(transpose(sol))
names!(df,[Symbol(reac) for reac in reactants])
CSV.write("data/results_gas_jac.csv",df)
df
#@profile run_simulation()
#open("prof.txt", "w") do s
#    Profile.print(IOContext(s, :displaysize => (1000, 500)))
#end
#using ProfileView
#ProfileView.view()
#Profile.clear()
#using Plots
#plot(log10.(sol[reactants2ind["BUT1ENE"],:]))
#plot!(log10.(sol[reactants2ind["C2H5O2"],:]))