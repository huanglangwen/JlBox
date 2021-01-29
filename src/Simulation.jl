function make_odefun(::GasConfig, solverconfig::SolverConfig, len_y::Int, jac_prototype)
    if solverconfig.diff_method == "gas_v"
        u = zeros(len_y)
        return DiffEqBase.ODEFunction(dydt!; jac_prototype = DiffEqOperators.AnalyticalJacVecOperator(gas_jac_v!, u))
    else
        return DiffEqBase.ODEFunction(dydt!; jac = gas_jac!, jac_prototype = jac_prototype)
    end
end

function make_odefun(::AerosolConfig, solverconfig::SolverConfig, len_y::Int, jac_prototype)
    if solverconfig.diff_method == "fine_seeding_v"
        u = zeros(len_y)
        return DiffEqBase.ODEFunction(dydt_aerosol!; jac_prototype = DiffEqOperators.AnalyticalJacVecOperator(aerosol_jac_v_fine_seeding!, u))
    else
        return DiffEqBase.ODEFunction(dydt_aerosol!; jac = select_jacobian(solverconfig.diff_method, len_y), jac_prototype = jac_prototype)
    end
end

function run_simulation(config::JlBoxConfig, solverconfig::SolverConfig)
    param_dict, reactants2ind, y_init = prepare(config, solverconfig)
    showconfig(config)
    showconfig(solverconfig, config.io)
    len_y = length(y_init)
    odefun = make_odefun(config, solverconfig, len_y, param_dict["sparsity"])
    odeprob = DiffEqBase.ODEProblem{true}(odefun,y_init,(0.,config.simulation_time),param_dict)
    @time sol = DiffEqBase.solve(odeprob, solverconfig.solver, reltol = solverconfig.reltol, abstol = solverconfig.abstol,
        tstops = 0:config.batch_step:config.simulation_time, saveat = config.batch_step, # save_everystep=true,
        dt = solverconfig.dtinit,#1.0e-6, #Initial step-size
        dtmax= solverconfig.dtmax,
        max_order = 5,
        max_convergence_failures = 1000,
        callback = solverconfig.positiveness ? DiffEqCallbacks.PositiveDomain(y_init) : nothing
        #isoutofdomain=(u,p,t) -> any(x -> x < 0, u)
        #progress=true
    )
    return sol, reactants2ind, param_dict
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

function run_simulation_aerosol_adjoint(aerosolconfig::AerosolConfig,aerosolsolverconfig::SolverConfig,adjointconfig::AdjointConfig)
    io = adjointconfig.io
    if isfile(joinpath(@__DIR__,"../data/aerosol_sol.store")) && adjointconfig.use_cache
        println(io, "Found cached aerosol simulation result")
        param_dict,_,y_init=prepare(aerosolconfig, aerosolsolverconfig)
        dy_dt_gas_matrix=zeros(Real,(param_dict["num_reactants"],aerosolconfig.num_bins))
        param_dict["dy_dt_gas_matrix"]=dy_dt_gas_matrix
        make_odefun(aerosolconfig, aerosolsolverconfig, length(y_init), param_dict["sparsity"])
        println(io, "Loading cache")
        sol=deserialize(joinpath(@__DIR__,"../data/aerosol_sol.store"))
    else
        println(io, "No caching, start aerosol simulation")
        sol,_,param_dict=run_simulation(aerosolconfig, aerosolsolverconfig)
        println(io, "Caching solution")
        serialize(joinpath(@__DIR__,"../data/aerosol_sol.store"),sol)
    end
    num_reactants,num_reactants_condensed,num_eqns=[param_dict[i] for i in ["num_reactants","num_reactants_condensed","num_eqns"]]
    println(io, "Preparing Adjoint Problem")
    t0,tF=(0.,aerosolconfig.simulation_time)
    tspan_adj=(tF,t0)
    len_y=num_reactants+aerosolconfig.num_bins*num_reactants_condensed
    mw_array=param_dict["y_mw"]
    lambda_init=zeros(Float64,(1,len_y))#DiffEq.jl version seems incorrect
    SOA_mass_jac!(lambda_init,mw_array,num_reactants,num_reactants_condensed,aerosolconfig.num_bins)#adopting KPP paper I
    param_dict["sol"]=sol
    param_dict["jac_mtx"]= aerosolsolverconfig.sparse ? spzeros(len_y, len_y) : zeros(Float64,(len_y,len_y))
    param_dict["Current_iter"]=0
    param_dict["ShowIterPeriod"]=5
    param_dict["Simulation_type"]="adjoint"
    param_dict["jacobian!"]=select_jacobian(adjointconfig.diff_method,len_y)
    param_dict["adjointconfig"] = adjointconfig
    odefun_adj=ODEFunction(sensitivity_adjoint_dldt!,jac=sensitivity_adjoint_jac!)
    prob_adj=ODEProblem{true}(odefun_adj,reshape(lambda_init, : ),tspan_adj,param_dict)
    println(io, "Using solver: $(typeof(adjointconfig.adjoint_solver))")
    println(io, "Using Jacobian: $(adjointconfig.diff_method)")
    println(io, "Reltol: $(adjointconfig.reltol), Abstol: $(adjointconfig.abstol)")
    @debug "Solving Adjoint Problem"
    lambda_sol=solve(prob_adj,adjointconfig.adjoint_solver,reltol=adjointconfig.reltol,abstol=adjointconfig.abstol,#Rodas5(autodiff=false)
                     tstops=aerosolconfig.simulation_time:-aerosolconfig.batch_step:0.,saveat=-aerosolconfig.batch_step,
                     dt=-1e-6,dtmax=50.0,max_order=5,max_convergence_failures=1000)
    @debug "Preparing Integration"
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
    @debug "Strating Integration"
    for i in 1:num_tstops-1
        println(io, "Integrating from $(tstops[i]) to $(tstops[i+1])")
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
    SOA_mass_jac!(dSOA_dy,mw_array,num_reactants,num_reactants_condensed,config.num_bins)
    return reshape(dSOA_dy * reshape(S,(y_len,num_eqns)),num_eqns)
end