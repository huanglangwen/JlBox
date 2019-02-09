include("JlBoxModule.jl")
using .Compute:run_simulation

#Profile.init(n = 10^7, delay = 5.)
@time sol,reactants2ind=run_simulation()

#@profile run_simulation()
#open("prof.txt", "w") do s
#    Profile.print(IOContext(s, :displaysize => (1000, 500)))
#end
#using ProfileView
#ProfileView.view()
#Profile.clear()
using Plots
plot(log10.(sol[reactants2ind["BUT1ENE"],:]))
plot!(log10.(sol[reactants2ind["C2H5O2"],:]))