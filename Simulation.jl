include("JlBoxModule.jl")
using .Compute:generate_codes,loss_gain!,tspan,batch_step,simulation_time,start_time,H2O,temp
using DifferentialEquations
using ForwardDiff
using Printf

println("Generating Codes")
reactants2ind,reactants_initial,dy,rate_values,J,stoich_mtx,stoich_list,RO2_inds,num_eqns,num_reactants=generate_codes()
include("generated_code.jl")
println("Solving ODE")

function dydt_jac!(J,reactants,p,t)
    ForwardDiff.jacobian!(J,x->dydt!(x,p,t),reactants)
    nothing
end

dydtf=ODEFunction(dydt!;jac=dydt_jac!)
prob = ODEProblem{false}(dydtf,reactants_initial,tspan,
                         (dy,rate_values,J,stoich_mtx,stoich_list,RO2_inds,num_eqns,num_reactants)
                        #(dy,rate_values,rate_prods,J,RO2_inds,num_eqns,num_reactants))
                        #(dy,rate_values,J,stoich_mtx,stoich_list,RO2_inds,num_eqns,num_reactants)
                        )
sol = solve(prob,CVODE_BDF(linear_solver=:Dense),reltol=1e-6,abstol=1.0e-3,#
            tstops=0:batch_step:simulation_time,saveat=batch_step,# save_everystep=true,
            dt=1.0e-6, #Initial step-size
            dtmax=100.0,
            max_order = 5,
            max_convergence_failures = 1000,#1000
            #progress=false #Juno Progressbar
            )

#using Profile
#Profile.init(n = 10^7, delay = 5.)
#@time sol,reactants2ind=run_simulation()

#@profile run_simulation()
#open("prof.txt", "w") do s
#    Profile.print(IOContext(s, :displaysize => (1000, 500)))
#end
#using ProfileView
#ProfileView.view()
#Profile.clear()
#using Plots
plot(log10.(sol[reactants2ind["BUT1ENE"],:]))
plot!(log10.(sol[reactants2ind["C2H5O2"],:]))