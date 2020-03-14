function run_simulation_gas(config)
    #read_configure!("Configure_gas.jl")
    param_dict,reactants2ind=prepare_gas(config)
    num_reactants=param_dict["num_reactants"]
    reactants_initial=zeros(Float64,num_reactants)
    for (k,v) in config.reactants_initial_dict
        reactants_initial[reactants2ind[k]]=v*config.Cfactor#pbb to molcules/cc
    end
    print("Using solver: ");println(typeof(config.solver))
    @printf("Reltol: %.3e, Abstol: %.3e\n",config.reltol,config.abstol)
    println("Solving ODE")
    param_dict["Current_iter"]=0
    param_dict["ShowIterPeriod"]=100
    param_dict["Simulation_type"]="gas"
    #odefun=ODEFunction(dydt!; jac=gas_jac!)
    if config.use_jacobian
        odefun=ODEFunction(dydt!; jac=gas_jac!)
        prob = ODEProblem{true}(odefun,reactants_initial,config.tspan,param_dict)
        param_dict["ShowIterPeriod"]=5
    else
        prob = ODEProblem{true}(dydt!,reactants_initial,config.tspan,param_dict)
        param_dict["ShowIterPeriod"]=500
    end
    @time sol = solve(prob,config.solver,reltol=config.reltol,abstol=config.abstol,
                tstops=0:config.batch_step:config.simulation_time,saveat=config.batch_step,# save_everystep=true,
                dt=1.0e-6, #Initial step-size
                dtmax=100.0,
                max_order = 5,
                max_convergence_failures = 1000,
                #callback=PositiveDomain(reactants_initial,abstol=1.0e-3)
                callback=config.positiveness ? PositiveDomain(reactants_initial) : nothing
                #isoutofdomain=(u,p,t) -> any(x -> x < 0, u)
                #progress=true
                )
    return sol,reactants2ind
end

function run_simulation_gas_sparse(solver, config, jac_prototype, param_dict, reactants2ind)
    num_reactants=param_dict["num_reactants"]
    reactants_initial=zeros(Float64,num_reactants)
    for (k,v) in config.reactants_initial_dict
        reactants_initial[reactants2ind[k]]=v*config.Cfactor#pbb to molcules/cc
    end
    print("Using solver: ");println(typeof(solver))
    @printf("Reltol: %.3e, Abstol: %.3e\n",config.reltol,config.abstol)
    println("Solving ODE")
    param_dict["Current_iter"]=0
    param_dict["Simulation_type"]="gas"
    param_dict["ShowIterPeriod"]=10
    odefun=ODEFunction(dydt!; jac=gas_jac!, jac_prototype=jac_prototype)
    #odefun=ODEFunction(dydt!; jac_prototype=DiffEqOperators.AnalyticalJacVecOperator(gas_jac_v!,reactants_initial,param_dict,0.0))
    prob = ODEProblem{true}(odefun,reactants_initial,config.tspan,param_dict)
    @time sol = solve(prob,solver,reltol=config.reltol,abstol=config.abstol,
                tstops=0:config.batch_step:config.simulation_time,saveat=config.batch_step,# save_everystep=true,
                dt=1.0e-6, #Initial step-size
                dtmax=100.0,
                max_order = 5,
                max_convergence_failures = 1000,
                #callback=PositiveDomain(reactants_initial,abstol=1.0e-3)
                callback=config.positiveness ? PositiveDomain(reactants_initial,abstol=1e2,scalefactor=0.9) : nothing
                #isoutofdomain=(u,p,t) -> any(x -> x < 0, u)
                #progress=true
                )
    return sol,reactants2ind
end

function select_jacobian(diff_method, len_y)
    if diff_method == "finite"
        jac_cache = FiniteDiff.JacobianCache(zeros(Float64,len_y),zeros(Float64,len_y),Val{:forward},Float64,Val{true})
        jac! = (jac_mtx,y,p,t)->FiniteDiff.finite_difference_jacobian!(jac_mtx,(dydt,y)->dydt_aerosol!(dydt,y,p,t),y,jac_cache)
    elseif diff_method == "coarse_seeding"
        jac! = aerosol_jac_coarse_seeding!
    elseif diff_method == "fine_seeding"
        jac! = aerosol_jac_fine_seeding!
    elseif diff_method == "analytical"
        jac! = aerosol_jac_fine_analytical!
    else
        println("ERROR: can't recognize diff type: ",diff_method)
        exit(-1)
    end
    jac!
end

function run_simulation_aerosol(config)
    #read_configure!("Configure_aerosol.jl")
    param_dict,reactants2ind,y_cond=prepare_aerosol(config)
    num_reactants,num_reactants_condensed=[param_dict[i] for i in ["num_reactants","num_reactants_condensed"]]
    dy_dt_gas_matrix=zeros(Real,(num_reactants,config.num_bins))
    #dy_dt=zeros(Real,num_reactants+num_reactants_condensed*num_bins)
    param_dict["dy_dt_gas_matrix"]=dy_dt_gas_matrix
    #param_dict["dydt"]=dy_dt
    param_dict["Current_iter"]=0
    param_dict["Simulation_type"]="aerosol"
    len_y=num_reactants+num_reactants_condensed*config.num_bins
    y_init=zeros(Float64,len_y)
    for (k,v) in config.reactants_initial_dict
        y_init[reactants2ind[k]]=v*config.Cfactor#pbb to molcules/cc
    end
    y_init[num_reactants+1:num_reactants+config.num_bins*num_reactants_condensed]=y_cond[1:config.num_bins*num_reactants_condensed]
    print("Using solver: ");println(typeof(config.solver))
    print("Using Jacobian: ");println(config.diff_method)
    @printf("Reltol: %.3e, Abstol: %.3e\n",config.reltol,config.abstol)
    println("Solving ODE")
    if config.use_jacobian
        odefun=ODEFunction(dydt_aerosol!; jac=select_jacobian(config.diff_method,len_y))
        prob = ODEProblem{true}(odefun,y_init,config.tspan,param_dict)
        param_dict["ShowIterPeriod"]=5
    else
        prob = ODEProblem{true}(dydt_aerosol!,y_init,config.tspan,param_dict)
        param_dict["ShowIterPeriod"]=500
    end
    sol = solve(prob,config.solver,reltol=config.reltol,abstol=config.abstol,
                tstops=0:config.batch_step:config.simulation_time,saveat=config.batch_step,# save_everystep=true,
                dt=1.0e-6, #Initial step-size
                dtmax=100.0,
                max_order = 5,
                max_convergence_failures = 1000,
                callback=config.positiveness ? PositiveDomain(y_init) : nothing
                #isoutofdomain=(u,p,t) -> any(x -> x < 0, u)
                )
    sol_mtx=transpose(sol)
    aerosol_mtx=sol_mtx[1:end,num_reactants+1:num_reactants+config.num_bins*num_reactants_condensed]
    t_length=size(aerosol_mtx)[1]
    mw_array=param_dict["y_mw"]
    SOA_array=[sum((sum(reshape(aerosol_mtx[i,1:end],(num_reactants_condensed,config.num_bins))
                               ,dims=2).*mw_array./config.NA)[1:end-1]#exclude H2O at the end
                  ) for i in 1:t_length]*1E12

    return sol,reactants2ind,SOA_array,num_reactants,param_dict
end

function rates_from_sol(p::Dict,t::Real)
    rate_values,J,RO2_inds,num_eqns,num_reactants=[p[ind] for ind in ["rate_values","J","RO2_inds","num_eqns","num_reactants"]]
    evaluate_rates_fun=p["evaluate_rates!"]
    config=p["config"]
    time_of_day_seconds=config.start_time+t
    sol=p["sol"]
    reactants=sol(t)[1:num_reactants]
    RO2=sum(reactants[RO2_inds])
    Base.invokelatest(evaluate_rates_fun,time_of_day_seconds,RO2,config.H2O,config.temp,rate_values,J)# =>ratevalues
    rate_values
end

function run_simulation_aerosol_adjoint(aerosolconfig,adjointconfig)
    #read_configure!("Configure_aerosol.jl")
    if isfile("../data/aerosol_sol.store") && adjointconfig.use_cache
        println("Found cached aerosol simulation result")
        #read_configure!("Configure_aerosol.jl")
        param_dict,_,_=prepare_aerosol(aerosolconfig)
        dy_dt_gas_matrix=zeros(Real,(param_dict["num_reactants"],aerosolconfig.num_bins))
        param_dict["dy_dt_gas_matrix"]=dy_dt_gas_matrix
        odefun=ODEFunction(dydt_aerosol!; jac=aerosol_jac_fine_seeding!)
        println("Loading cache")
        sol=deserialize("../data/aerosol_sol.store")
    else
        println("No caching, start aerosol simulation")
        sol,_,_,_,param_dict=run_simulation_aerosol(aerosolconfig)
        println("Caching solution")
        serialize("../data/aerosol_sol.store",sol)
    end
    num_reactants,num_reactants_condensed,num_eqns=[param_dict[i] for i in ["num_reactants","num_reactants_condensed","num_eqns"]]
    println("Preparing Adjoint Problem")
    t0,tF=aerosolconfig.tspan
    tspan_adj=(tF,t0)
    len_y=num_reactants+aerosolconfig.num_bins*num_reactants_condensed
    mw_array=param_dict["y_mw"]
    lambda_init=zeros(Float64,(1,len_y))#DiffEq.jl version seems incorrect
    SOA_mass_jac!(lambda_init,mw_array,aerosolconfig.NA,num_reactants,num_reactants_condensed,aerosolconfig.num_bins)#adopting KPP paper I
    #println(lambda_init)
    param_dict["sol"]=sol
    param_dict["jac_mtx"]=zeros(Float64,(len_y,len_y))
    param_dict["Current_iter"]=0
    param_dict["ShowIterPeriod"]=5
    param_dict["Simulation_type"]="adjoint"
    param_dict["jacobian!"]=select_jacobian(adjointconfig.diff_method,len_y)
    odefun_adj=ODEFunction(sensitivity_adjoint_dldt!,jac=sensitivity_adjoint_jac!)
    prob_adj=ODEProblem{true}(odefun_adj,reshape(lambda_init, : ),tspan_adj,param_dict)
    print("Using solver: ");println(typeof(adjointconfig.adjoint_solver))
    print("Using Jacobian: ");println(adjointconfig.diff_method)
    @printf("Reltol: %.3e, Abstol: %.3e\n",adjointconfig.reltol,adjointconfig.abstol)
    println("Solving Adjoint Problem")
    lambda_sol=solve(prob_adj,adjointconfig.adjoint_solver,reltol=adjointconfig.reltol,abstol=adjointconfig.abstol,#Rodas5(autodiff=false)
                     tstops=aerosolconfig.simulation_time:-aerosolconfig.batch_step:0.,saveat=-aerosolconfig.batch_step,
                     dt=-1e-6,dtmax=50.0,max_order=5,max_convergence_failures=1000)
    println("Preparing Integration")
    tstops=[t for t in 0:aerosolconfig.batch_step:aerosolconfig.simulation_time]
    num_tstops=length(tstops)
    stoich_mtx=param_dict["stoich_mtx"]
    stoich_list=param_dict["stoich_list"]
    reactants_list=param_dict["reactants_list"]
    dSOA_mass_drate=zeros(Float64,(num_eqns,num_tstops))
    dSOA_mass_percentk=zeros(Float64,(num_eqns,num_tstops))
    dgpdt=function (t)
        y_gas=sol(t)[1:num_reactants]
        loss_gain_drate_mtx=zeros(Float64,(num_reactants,num_eqns))
        loss_gain_drate_values!(num_reactants,num_eqns,y_gas,stoich_mtx,stoich_list,reactants_list,loss_gain_drate_mtx)
        lambda=lambda_sol(t)[1:num_reactants]
        return lambda' * loss_gain_drate_mtx
    end 
    println("Strating Integration")
    for i in 1:num_tstops-1
        @printf("Integrating from %.0f to %.0f\n",tstops[i],tstops[i+1])
        dSOA_mass_drate[1:num_eqns,i+1]=dSOA_mass_drate[1:num_eqns,i]+reshape(quadgk(dgpdt,tstops[i],tstops[i+1])[1],(num_eqns,1))#quadgk->(val,err) ignore error value
        rate_values=rates_from_sol(param_dict,tstops[i+1])
        dSOA_mass_percentk[1:num_eqns,i+1]=0.01*rate_values.*dSOA_mass_drate[1:num_eqns,i+1]
    end
    return dSOA_mass_drate,dSOA_mass_percentk
end

function sensitivity_mtx2dSOA(S,t::Real,integrator)
    p=integrator.p
    config=p["config"]
    mw_array,num_reactants,num_reactants_condensed,num_eqns=[p[i] for i in ["y_mw","num_reactants","num_reactants_condensed","num_eqns"]]
    y_len=num_reactants+config.num_bins*num_reactants_condensed
    dSOA_dy=zeros(Float64,(1,y_len))
    SOA_mass_jac!(dSOA_dy,mw_array,config.NA,num_reactants,num_reactants_condensed,config.num_bins)
    #println(dSOA_dy)
    println(size(S))
    println(S[1:100])
    return reshape(dSOA_dy * reshape(S,(y_len,num_eqns)),num_eqns)
end