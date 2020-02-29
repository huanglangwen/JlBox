"""
loss_gainl!(num_reactants,num_eqns,reactants,stoich_mtx,
            stoich_list,reactants_list,rate_values,dy_dt)

This function is intended to compute the value of RHS function
for gas kinetics and store them in `dy_dt`. `stoich_list`,
`reactants_list` are both produced from `stoich_mtx` by
`mk_reactants_list(...)` to reduce access overhead of CSC
matrices. 
"""
function loss_gain!(num_reactants::Int,num_eqns::Int,
                   reactants::Array{<:Real,1},#num_reactants
                   stoich_mtx::SparseMatrixCSC{Float64,Int64},#num_reactants*num_eqns
                   stoich_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, both reac and prod
                   reactants_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, only reac
                   rate_values::Array{<:Real,1},#num_eqns
                   dydt::Array{<:Real,1}#num_reactants
                   )
    lossgain_mtx=similar(stoich_mtx)#spzeros(num_reactants,num_eqns)
    @inbounds for eqn_ind in 1:num_eqns
        prod=rate_values[eqn_ind]
        num_reacs,stoichvec,indvec=reactants_list[eqn_ind]
        num_stoichs,_,stoich_indvec=stoich_list[eqn_ind]
        @inbounds for i in 1:num_reacs
            reactant_ind=indvec[i]
            stoich=stoichvec[i]
            #stoich=stoich_mtx[reactant_ind,eqn_ind]
            prod*=reactants[reactant_ind]^stoich#reactants_list come from reactants_mtx (for catalyse A+B=A+C)
            #if stoich>0
            #    prod*=reactants[reactant_ind]^stoich
            #end
        end
        #lossgain_mtx[:,eqn_ind].=stoich_mtx[:,eqn_ind].*prod #!!!!! HUGE PERFORMANCE COST
        #for reactant_ind in reactant_inds
        @inbounds for i in 1:num_stoichs
            reactant_ind=stoich_indvec[i]
            lossgain_mtx[reactant_ind,eqn_ind]=stoich_mtx[reactant_ind,eqn_ind]*prod
        end
    end
    fill!(dydt,0.)
    @inbounds for i = 1:length(lossgain_mtx.rowval)
        reactant_ind=lossgain_mtx.rowval[i]
        val=lossgain_mtx.nzval[i]
        dydt[reactant_ind]-=val
    end
    #is,js,vs=findnz(lossgain_mtx)
    #lossgain_mtx_T=sparse(js,is,vs,num_eqns,num_reactants)#num_eqns*num_reactants
    #for reactant_ind in 1:num_reactants
    #    dydt[reactant_ind]=sum(nonzeros(lossgain_mtx_T[:,reactant_ind]))*(-1)#dydt negative for reactants, positive for products 
    #end #*reactants[reactant_ind]=>wrong!!!
    return dydt
end

function dydt!(dydt,reactants::Array{<:Real,1},p::Dict,t::Real)
    #dy,rate_values,J,stoich_mtx,stoich_list,reactants_list,RO2_inds,num_eqns,num_reactants=p
    rate_values,J,stoich_mtx,stoich_list,reactants_list,RO2_inds,num_eqns,num_reactants=
        [p[ind] for ind in 
            ["rate_values","J","stoich_mtx","stoich_list","reactants_list","RO2_inds",
             "num_eqns","num_reactants"]
        ]
    config=p["config"]
    evaluate_rates_fun=p["evaluate_rates!"]
    time_of_day_seconds=config.start_time+t
    RO2=sum(reactants[RO2_inds])
    Base.invokelatest(evaluate_rates_fun,time_of_day_seconds,RO2,config.H2O,config.temp,rate_values,J)# =>ratevalues
    loss_gain!(num_reactants,num_eqns,reactants,stoich_mtx,stoich_list,reactants_list,rate_values,dydt)
    #loss_gain_static!(num_reactants,num_eqns,reactants,rate_values,rate_prods,dy)
    if p["Simulation_type"]=="gas"
        p["Current_iter"]+=1
        citer=p["Current_iter"]
        if citer%(p["ShowIterPeriod"])==0
            @printf("Current Iteration: %d, time_step: %e\n",citer,t)
        end
    end
    nothing#return dydt
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

function dydt_aerosol!(dy_dt,y::Array{<:Real,1},p::Dict,t::Real)
    num_reactants,num_reactants_condensed=[p[i] for i in ["num_reactants","num_reactants_condensed"]]
    include_inds,dy_dt_gas_matrix,N_perbin=[p[i] for i in ["include_inds","dy_dt_gas_matrix","N_perbin"]]
    config=p["config"]
    mw_array,density_array,gamma_gas,alpha_d_org,DStar_org,Psat=[p[i] for i in ["y_mw","y_density_array","gamma_gas","alpha_d_org","DStar_org","Psat"]]
    y_core,core_mass_array=[p[i] for i in ["y_core","core_mass_array"]]
    y_gas=y[1:num_reactants]#view(xs,lo:hi) passes ref instead of copy
    dydt!(dy_dt,y_gas,p,t)
    C_g_i_t=y[include_inds]
    _,total_SOA_mass=Partition!(y,dy_dt,dy_dt_gas_matrix,C_g_i_t,
        config.num_bins,num_reactants,num_reactants_condensed,include_inds,
        mw_array,density_array,gamma_gas,alpha_d_org,DStar_org,Psat,N_perbin,
        config.core_dissociation,y_core,core_mass_array,config.core_density_array,
        config.NA,config.sigma,config.R_gas,config.temp)
    if p["Simulation_type"]=="aerosol"
        p["Current_iter"]+=1
        citer=p["Current_iter"]
        if citer%(p["ShowIterPeriod"])==0
            @printf("Current Iteration: %d, time_step: %e, SOA(ug/m3): %e\n",citer,t,total_SOA_mass)
            #println("Sum(dy_dt[num_reacs+1:end])=",sum(dy_dt[num_reactants+1:end]))
            #println("Sum(y[num_reacs+1:end])=",sum(y[num_reactants+1:end]))
        end
    end
    nothing#return dy_dt
end

function aerosol_jac!(jac_mtx,y::Array{Float64,1},p::Dict,t::Real)
    num_reactants,num_reactants_condensed=[p[i] for i in ["num_reactants","num_reactants_condensed"]]
    rate_values,J,RO2_inds=[p[i] for i in ["rate_values","J","RO2_inds"]]
    evaluate_rates_fun=p["evaluate_rates!"]
    config=p["config"]
    time_of_day_seconds=config.start_time+t
    RO2=sum(y[RO2_inds])
    Base.invokelatest(evaluate_rates_fun,time_of_day_seconds,RO2,config.H2O,config.temp,rate_values,J)
    gas_jac!(jac_mtx,y,p,t)
    include_inds,dy_dt_gas_matrix,N_perbin=[p[i] for i in ["include_inds","dy_dt_gas_matrix","N_perbin"]]
    mw_array,density_array,gamma_gas,alpha_d_org,DStar_org,Psat=[p[i] for i in ["y_mw","y_density_array","gamma_gas","alpha_d_org","DStar_org","Psat"]]
    y_core,core_mass_array=[p[i] for i in ["y_core","core_mass_array"]]
    C_g_i_t=y[include_inds]
    Partition_jac_AD!(jac_mtx,y,C_g_i_t,
        config.num_bins,num_reactants,num_reactants_condensed,include_inds,
        mw_array,density_array,gamma_gas,alpha_d_org,DStar_org,Psat,N_perbin,
        config.core_dissociation,y_core,core_mass_array,config.core_density_array,
        config.NA,config.sigma,config.R_gas,config.temp)
    nothing
end

function aerosol_jac_seeding!(jac_mtx,y::Array{Float64,1},p::Dict,t::Real)
    num_reactants,num_reactants_condensed=[p[i] for i in ["num_reactants","num_reactants_condensed"]]
    rate_values,J,RO2_inds=[p[i] for i in ["rate_values","J","RO2_inds"]]
    evaluate_rates_fun=p["evaluate_rates!"]
    config=p["config"]
    time_of_day_seconds=config.start_time+t
    RO2=sum(y[RO2_inds])
    Base.invokelatest(evaluate_rates_fun,time_of_day_seconds,RO2,config.H2O,config.temp,rate_values,J)
    
    include_inds,dy_dt_gas_matrix,N_perbin=[p[i] for i in ["include_inds","dy_dt_gas_matrix","N_perbin"]]
    mw_array,density_array,gamma_gas,alpha_d_org,DStar_org,Psat=[p[i] for i in ["y_mw","y_density_array","gamma_gas","alpha_d_org","DStar_org","Psat"]]
    y_core,core_mass_array=[p[i] for i in ["y_core","core_mass_array"]]
    
    partition_dydt_fun=function (dy_dt,y)
        C_g_i_t=y[include_inds]
        Partition!(y,dy_dt,dy_dt_gas_matrix,C_g_i_t,
        config.num_bins,num_reactants,num_reactants_condensed,include_inds,
        mw_array,density_array,gamma_gas,alpha_d_org,DStar_org,Psat,N_perbin,
        config.core_dissociation,y_core,core_mass_array,config.core_density_array,
        config.NA,config.sigma,config.R_gas,config.temp)
    end
    dy_dt=zeros(Real,length(y))
    ForwardDiff.jacobian!(jac_mtx,partition_dydt_fun, dy_dt, y)
    gas_jac!(jac_mtx,y,p,t)
    nothing
end

function sensitivity_adjoint_jac!(jac_mtx,lambda,p,t)
    jacobian_from_sol!(p,t)#jacobian_from_sol!(p,t)
    jac_mtx.=(-1).*transpose(p["jac_mtx"])#IMPORTANT jacobian should be the transpose of the original one 
    # since dldt=g(t)-l*J, for ith element in l and jth element in dldt appears at ith line and jth col in the Jacobian matrix
    nothing
end

function jacobian_from_sol!(p::Dict,t::Real)
    sol=p["sol"]
    y=sol(t)
    jac_mtx=p["jac_mtx"]
    fill!(jac_mtx,0.)
    diff=p["Diff_method"]
    if diff=="finite"
        jac_cache=p["jac_cache"]
        FiniteDiff.finite_difference_jacobian!(jac_mtx,(dydt,y)->dydt_aerosol!(dydt,y,p,t),y,jac_cache)
    elseif diff=="dual"
        aerosol_jac_seeding!(jac_mtx,y,p,t)
    elseif diff=="analytical"
        aerosol_jac!(jac_mtx,y,p,t)
    else
        println("WARNING: can't recognize diff type: ",diff)
    end
    nothing
end

function sensitivity_adjoint_dldt!(dldt,lambda,p,t)
    jacobian_from_sol!(p,t)#jacobian_from_sol!(p,t)
    jac_mtx=p["jac_mtx"]
    dldt.= reshape(- lambda' * jac_mtx, :)#adopting KPP paper I
    p["Current_iter"]+=1
    citer=p["Current_iter"]
    if citer%(p["ShowIterPeriod"])==0
        num_reactants=p["num_reactants"]
        @printf("Current Iteration: %d, time_step: %e, sum(lambda_gas): %e, sum(dldt_gas): %e, sum(lambda): %e\n",citer,t,sum(lambda[1:num_reactants]),sum(dldt[1:num_reactants]),sum(lambda))
        #println(sum(jac_mtx[:,1:num_reactants],dims=1))
    end
    nothing
end

function run_simulation_aerosol(config;use_jacobian::Bool)
    #read_configure!("Configure_aerosol.jl")
    param_dict,reactants2ind,y_cond=prepare_aerosol(config)
    num_reactants,num_reactants_condensed=[param_dict[i] for i in ["num_reactants","num_reactants_condensed"]]
    dy_dt_gas_matrix=zeros(Real,(num_reactants,config.num_bins))
    #dy_dt=zeros(Real,num_reactants+num_reactants_condensed*num_bins)
    param_dict["dy_dt_gas_matrix"]=dy_dt_gas_matrix
    #param_dict["dydt"]=dy_dt
    param_dict["Current_iter"]=0
    param_dict["Simulation_type"]="aerosol"
    y_init=zeros(Float64,num_reactants+num_reactants_condensed*config.num_bins)
    for (k,v) in config.reactants_initial_dict
        y_init[reactants2ind[k]]=v*config.Cfactor#pbb to molcules/cc
    end
    y_init[num_reactants+1:num_reactants+config.num_bins*num_reactants_condensed]=y_cond[1:config.num_bins*num_reactants_condensed]
    println("Solving ODE")
    if use_jacobian
        odefun=ODEFunction(dydt_aerosol!; jac=aerosol_jac!)
        prob = ODEProblem{true}(odefun,y_init,config.tspan,param_dict)
        param_dict["ShowIterPeriod"]=5
    else
        prob = ODEProblem{true}(dydt_aerosol!,y_init,config.tspan,param_dict)
        param_dict["ShowIterPeriod"]=500
    end
    sol = solve(prob,config.solver,reltol=1e-4,abstol=1.0e-2,
                tstops=0:config.batch_step:config.simulation_time,saveat=config.batch_step,# save_everystep=true,
                dt=1.0e-6, #Initial step-size
                dtmax=100.0,
                max_order = 5,
                max_convergence_failures = 1000,
                #callback=PositiveDomain(y_init,abstol=1.0e-2)
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

function run_simulation_aerosol_adjoint(config)
    #read_configure!("Configure_aerosol.jl")
    if isfile("../data/aerosol_sol.store")
        println("Found cached aerosol simulation result")
        #read_configure!("Configure_aerosol.jl")
        param_dict,_,_=prepare_aerosol(config)
        dy_dt_gas_matrix=zeros(Real,(param_dict["num_reactants"],config.num_bins))
        param_dict["dy_dt_gas_matrix"]=dy_dt_gas_matrix
        odefun=ODEFunction(dydt_aerosol!; jac=aerosol_jac!)
        println("Loading cache")
        sol=deserialize("../data/aerosol_sol.store")
    else
        println("No caching, start aerosol simulation")
        sol,_,_,_,param_dict=run_simulation_aerosol(config,use_jacobian=true)
        println("Caching solution")
        serialize("../data/aerosol_sol.store",sol)
    end
    num_reactants,num_reactants_condensed,num_eqns=[param_dict[i] for i in ["num_reactants","num_reactants_condensed","num_eqns"]]
    println("Preparing Adjoint Problem")
    t0,tF=config.tspan
    tspan_adj=(tF,t0)
    len_y=num_reactants+config.num_bins*num_reactants_condensed
    mw_array=param_dict["y_mw"]
    lambda_init=zeros(Float64,(1,len_y))#DiffEq.jl version seems incorrect
    SOA_mass_jac!(lambda_init,mw_array,config.NA,num_reactants,num_reactants_condensed,config.num_bins)#adopting KPP paper I
    #println(lambda_init)
    param_dict["sol"]=sol
    param_dict["jac_mtx"]=zeros(Float64,(len_y,len_y))
    param_dict["Current_iter"]=0
    param_dict["ShowIterPeriod"]=5
    param_dict["Simulation_type"]="adjoint"
    param_dict["Diff_method"]="dual"
    param_dict["jac_cache"]=FiniteDiff.JacobianCache(zeros(Float64,len_y),Val{:forward},Float64,Val{true})
    odefun_adj=ODEFunction(sensitivity_adjoint_dldt!,jac=sensitivity_adjoint_jac!)
    prob_adj=ODEProblem{true}(odefun_adj,reshape(lambda_init, : ),tspan_adj,param_dict)
    println("Solving Adjoint Problem")
    lambda_sol=solve(prob_adj,config.solver,reltol=1e-8,abstol=1e-6,#Rodas5(autodiff=false)
                     tstops=config.simulation_time:-config.batch_step:0.,saveat=-config.batch_step,
                     dt=-1e-6,dtmax=50.0,max_order=5,max_convergence_failures=1000)
    println("Preparing Integration")
    tstops=[t for t in 0:config.batch_step:config.simulation_time]
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

function run_simulation_gas(config;use_jacobian::Bool=true)
    #read_configure!("Configure_gas.jl")
    param_dict,reactants2ind=prepare_gas(config)
    num_reactants=param_dict["num_reactants"]
    reactants_initial=zeros(Float64,num_reactants)
    for (k,v) in config.reactants_initial_dict
        reactants_initial[reactants2ind[k]]=v*config.Cfactor#pbb to molcules/cc
    end
    println("Solving ODE")
    param_dict["Current_iter"]=0
    param_dict["ShowIterPeriod"]=100
    param_dict["Simulation_type"]="gas"
    #odefun=ODEFunction(dydt!; jac=gas_jac!)
    if use_jacobian
        odefun=ODEFunction(dydt!; jac=gas_jac!)
        prob = ODEProblem{true}(odefun,reactants_initial,config.tspan,param_dict)
        param_dict["ShowIterPeriod"]=5
    else
        prob = ODEProblem{true}(dydt!,reactants_initial,config.tspan,param_dict)
        param_dict["ShowIterPeriod"]=500
    end
    @time sol = solve(prob,config.solver,reltol=1e-6,abstol=1.0e-3,
                tstops=0:config.batch_step:config.simulation_time,saveat=config.batch_step,# save_everystep=true,
                dt=1.0e-6, #Initial step-size
                dtmax=100.0,
                max_order = 5,
                max_convergence_failures = 1000,
                #callback=PositiveDomain(reactants_initial,abstol=1.0e-3)
                #isoutofdomain=(u,p,t) -> any(x -> x < 0, u)
                #progress=true
                )
    return sol,reactants2ind
end
