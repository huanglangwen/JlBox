#include("JlBoxModule.jl")
using ..Parse_eqn:parse_reactants,gen_evaluate_rates
using ..Optimize:constant_folding!,extract_constants!,generate_loss_gain,mk_reactants_list
using ..SizeDistributions:lognormal
using ..PropertyCalculation:Pure_component1,Pure_component2
using ..Partitioning:Partition!
using ..Jacobian:gas_jac!,Partition_jac!
using ..Sensitivity:SOA_mass_jac!,loss_gain_drate_values!
using DifferentialEquations
using DifferentialEquations:CVODE_BDF,CVODE_Adams
using DiffEqDiffTools
using ForwardDiff
using ForwardDiff:JacobianConfig
using OrdinaryDiffEq
using StaticArrays
using SparseArrays
using Printf
#using DiffEqSensitivity
using QuadGK
#using BSON: @save, @load
using Serialization

function loss_gain!(num_reactants::Int,num_eqns::Int,
                   reactants::Array{<:Real,1},#num_reactants
                   stoich_mtx::SparseMatrixCSC{Float64,Int64},#num_reactants*num_eqns
                   stoich_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, both reac and prod
                   reactants_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, only reac
                   rate_values::Array{<:Real,1},#num_eqns
                   dydt::Array{<:Real,1}#num_reactants
                   )
    lossgain_mtx=spzeros(num_reactants,num_eqns)
    for eqn_ind in 1:num_eqns
        prod=rate_values[eqn_ind]
        num_reacs,stoichvec,indvec=reactants_list[eqn_ind]
        num_stoichs,_,stoich_indvec=stoich_list[eqn_ind]
        for i in 1:num_reacs
            reactant_ind=indvec[i]
            stoich=stoichvec[i]
            #stoich=stoich_mtx[reactant_ind,eqn_ind]
            prod*=reactants[reactant_ind]^stoich#reactants_list come from reactants_mtx (for catalyse A+B=A+C)
            #if stoich>0
            #    prod*=reactants[reactant_ind]^stoich
            #end
        end
        #lossgain_mtx[:,eqn_ind]=stoich_mtx[:,eqn_ind]*prod #!!!!! HUGE PERFORMANCE COST
        #for reactant_ind in reactant_inds
        for i in 1:num_stoichs
            reactant_ind=stoich_indvec[i]
            lossgain_mtx[reactant_ind,eqn_ind]=stoich_mtx[reactant_ind,eqn_ind]*prod
        end
    end
    is,js,vs=findnz(lossgain_mtx)
    lossgain_mtx_T=sparse(js,is,vs,num_eqns,num_reactants)#num_eqns*num_reactants
    for reactant_ind in 1:num_reactants
        dydt[reactant_ind]=sum(nonzeros(lossgain_mtx_T[:,reactant_ind]))*(-1)#dydt negative for reactants, positive for products 
    end #*reactants[reactant_ind]=>wrong!!!
    return dydt
end

function dydt!(dydt,reactants::Array{<:Real,1},p::Dict,t::Real)
    #dy,rate_values,J,stoich_mtx,stoich_list,reactants_list,RO2_inds,num_eqns,num_reactants=p
    rate_values,J,stoich_mtx,stoich_list,reactants_list,RO2_inds,num_eqns,num_reactants=
        [p[ind] for ind in 
            ["rate_values","J","stoich_mtx","stoich_list","reactants_list","RO2_inds",
             "num_eqns","num_reactants"]
        ]
    evaluate_rates_fun=p["evaluate_rates!"]
    time_of_day_seconds=start_time+t
    RO2=sum(reactants[RO2_inds])
    Base.invokelatest(evaluate_rates_fun,time_of_day_seconds,RO2,H2O,temp,rate_values,J)# =>ratevalues
    loss_gain!(num_reactants,num_eqns,reactants,stoich_mtx,stoich_list,reactants_list,rate_values,dydt)
    #loss_gain_static!(num_reactants,num_eqns,reactants,rate_values,rate_prods,dy)
    #if p["Simulation_type"]=="gas"
    #    p["Current_iter"]+=1
    #    citer=p["Current_iter"]
    #    if citer%(p["ShowIterPeriod"])==0
    #        @printf("Current Iteration: %d, time_step: %e\n",citer,t)
    #    end
    #end
    nothing#return dydt
end

function dydt_aerosol!(dy_dt,y::Array{<:Real,1},p::Dict,t::Real)
    num_reactants,num_reactants_condensed=[p[i] for i in ["num_reactants","num_reactants_condensed"]]
    include_inds,dy_dt_gas_matrix,N_perbin=[p[i] for i in ["include_inds","dy_dt_gas_matrix","N_perbin"]]
    mw_array,density_array,gamma_gas,alpha_d_org,DStar_org,Psat=[p[i] for i in ["y_mw","y_density_array","gamma_gas","alpha_d_org","DStar_org","Psat"]]
    y_core,core_mass_array=[p[i] for i in ["y_core","core_mass_array"]]
    y_gas=y[1:num_reactants]#view(xs,lo:hi) passes ref instead of copy
    dydt!(dy_dt,y_gas,p,t)
    C_g_i_t=y[include_inds]
    _,total_SOA_mass=Partition!(y,dy_dt,dy_dt_gas_matrix,C_g_i_t,
        num_bins,num_reactants,num_reactants_condensed,include_inds,
        mw_array,density_array,gamma_gas,alpha_d_org,DStar_org,Psat,N_perbin,
        core_dissociation,y_core,core_mass_array,core_density_array,
        NA,sigma,R_gas,temp)
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
    time_of_day_seconds=start_time+t
    RO2=sum(y[RO2_inds])
    Base.invokelatest(evaluate_rates_fun,time_of_day_seconds,RO2,H2O,temp,rate_values,J)
    gas_jac!(jac_mtx,y,p,t)
    include_inds,dy_dt_gas_matrix,N_perbin=[p[i] for i in ["include_inds","dy_dt_gas_matrix","N_perbin"]]
    mw_array,density_array,gamma_gas,alpha_d_org,DStar_org,Psat=[p[i] for i in ["y_mw","y_density_array","gamma_gas","alpha_d_org","DStar_org","Psat"]]
    y_core,core_mass_array=[p[i] for i in ["y_core","core_mass_array"]]
    C_g_i_t=y[include_inds]
    Partition_jac!(jac_mtx,y,C_g_i_t,
        num_bins,num_reactants,num_reactants_condensed,include_inds,
        mw_array,density_array,gamma_gas,alpha_d_org,DStar_org,Psat,N_perbin,
        core_dissociation,y_core,core_mass_array,core_density_array,
        NA,sigma,R_gas,temp)
    nothing
end

function aerosol_jac_seeding!(jac_mtx,y::Array{Float64,1},p::Dict,t::Real)
    num_reactants,num_reactants_condensed=[p[i] for i in ["num_reactants","num_reactants_condensed"]]
    rate_values,J,RO2_inds=[p[i] for i in ["rate_values","J","RO2_inds"]]
    evaluate_rates_fun=p["evaluate_rates!"]
    time_of_day_seconds=start_time+t
    RO2=sum(y[RO2_inds])
    Base.invokelatest(evaluate_rates_fun,time_of_day_seconds,RO2,H2O,temp,rate_values,J)
    
    include_inds,dy_dt_gas_matrix,N_perbin=[p[i] for i in ["include_inds","dy_dt_gas_matrix","N_perbin"]]
    mw_array,density_array,gamma_gas,alpha_d_org,DStar_org,Psat=[p[i] for i in ["y_mw","y_density_array","gamma_gas","alpha_d_org","DStar_org","Psat"]]
    y_core,core_mass_array=[p[i] for i in ["y_core","core_mass_array"]]
    C_g_i_t=y[include_inds]
    partition_dydt_fun=function (dydt,y)
        Partition!(y,dy_dt,dy_dt_gas_matrix,C_g_i_t,
        num_bins,num_reactants,num_reactants_condensed,include_inds,
        mw_array,density_array,gamma_gas,alpha_d_org,DStar_org,Psat,N_perbin,
        core_dissociation,y_core,core_mass_array,core_density_array,
        NA,sigma,R_gas,temp)
    end
    dy_dt=zeros(Real,length(y))
    ForwardDiff.jacobian!(jac_mtx,partition_dydt_fun, dy_dt, y)
    gas_jac!(jac_mtx,y,p,t)
    nothing
end

function sensitivity_adjoint_jac!(jac_mtx,lambda,p,t)
    jacobian_from_sol_finitediff!(p,t)#jacobian_from_sol!(p,t)
    jac_mtx=(-1).*transpose(p["jac_mtx"])#IMPORTANT jacobian should be the transpose of the original one 
    # since dldt=g(t)-l*J, for ith element in l and jth element in dldt appears at ith line and jth col in the Jacobian matrix
    nothing
end

function jacobian_from_sol!(p::Dict,t::Real;finitediff=false)
    sol=p["sol"]
    y=sol(t)
    jac_mtx=p["jac_mtx"]
    fill!(jac_mtx,0.)
    if finitediff
        jac_cache=p["jac_cache"]
        DiffEqDiffTools.finite_difference_jacobian!(jac_mtx,(dydt,y)->dydt_aerosol!(dydt,y,p,t),y,jac_cache)
    else
        aerosol_jac!(jac_mtx,y,p,t)
    end
    nothing
end

function sensitivity_adjoint_dldt!(dldt,lambda,p,t)
    jacobian_from_sol!(p,t,finitediff=true)#jacobian_from_sol!(p,t)
    jac_mtx=p["jac_mtx"]
    #dSOA_dy=zeros(Float64,(1,num_reactants+num_bins*num_reactants_condensed))
    #SOA_mass_jac!(dSOA_dy,mw_array,NA,num_reactants,num_reactants_condensed,num_bins)

    dldt=-lambda*jac_mtx#adopting KPP paper I
    if isnan(sum(lambda))
        println(lambda)
    end
    p["Current_iter"]+=1
    citer=p["Current_iter"]
    if citer%(p["ShowIterPeriod"])==0
        num_reactants=p["num_reactants"]
        @printf("Current Iteration: %d, time_step: %e, sum(lambda_gas): %e, sum(dldt_gas): %e, sum(lambda): %e\n",citer,t,sum(lambda[1:num_reactants]),sum(dldt[1:num_reactants]),sum(lambda))
        #println(sum(jac_mtx[:,1:num_reactants],dims=1))
    end
    nothing
end

function prepare_gas()
    println("Parsing Reactants")
    stoich_mtx,reactants_mtx,RO2_inds,num_eqns,num_reactants,reactants2ind=parse_reactants(file)
    reactants_list=mk_reactants_list(num_reactants,num_eqns,reactants_mtx)
    stoich_list=mk_reactants_list(num_reactants,num_eqns,stoich_mtx)
    @printf("num_eqns: %d, num_reactants: %d\n",num_eqns,num_reactants)

    println("Generating evaluate_rates()")
    evaluate_rates_expr=gen_evaluate_rates(file)
    println("Done Generation")
    rate_values=zeros(Float64,num_eqns)
    rate_prods=zeros(Float64,num_eqns)
    J=zeros(Float64,62)
    #dydt=zeros(Float64,num_reactants)
    println("Performing constant folding")
    constant_folding!(evaluate_rates_expr,constantdict,rate_values);
    extract_constants!(evaluate_rates_expr);
    println("Evaluating evaluate_rates&loss_gain codes")
    eval(evaluate_rates_expr)
    param_dict=Dict("rate_values"=>rate_values,"J"=>J,"stoich_mtx"=>stoich_mtx,#"dydt"=>dydt,
                    "stoich_list"=>stoich_list,"reactants_list"=>reactants_list,"RO2_inds"=>RO2_inds,
                    "num_eqns"=>num_eqns,"num_reactants"=>num_reactants,"evaluate_rates!"=>evaluate_rates!)
    return param_dict,reactants2ind
end

function prepare_aerosol()
    param_dict,reactants2ind=prepare_gas()
    num_reactants=param_dict["num_reactants"]
    ind2reactants=Dict(reactants2ind[reac]=>reac for reac in keys(reactants2ind))
    species_names=[ind2reactants[ind] for ind=1:num_reactants]

    println("Calculating Partitioning Properties: Part1")
    pc1_dict=Pure_component1(num_reactants,species_names,vp_cutoff,temp,property_methods)
    
    println("Adding H2O")
    num_reactants+=1
    param_dict["num_reactants"]=num_reactants#not pc1_dict
    push!(pc1_dict["include_inds"],num_reactants)
    reactants2ind["H2O"]=num_reactants
    include_inds=pc1_dict["include_inds"]
    num_reactants_condensed=length(include_inds)
    sat_vap_water = exp(-0.58002206E4/temp+0.13914993E1-
        0.48640239E-1*temp+0.41764768E-4*(temp^2.0E0)-
        0.14452093E-7*(temp^3.0E0)+0.65459673E1*log(temp))#Pa
    push!(pc1_dict["y_density_array"],1000.0E0)#Append density of water to array [kg/m3]
    push!(pc1_dict["y_mw"],18.0E0)#Append mw of water to array [g/mol]
    push!(pc1_dict["Psat"],sat_vap_water*9.86923E-6)#Convert Pa to atm
    push!(pc1_dict["Delta_H"],40.66)
    Lv_water_vapour=2.5e3 # Latent heat of vapourisation of water [J/g] 
    push!(pc1_dict["Latent_heat_gas"],Lv_water_vapour)#Water vapour, taken from Paul Connolly's parcel model ACPIM

    println("Calculating Partitioning Properties: Part2")
    y_mw=pc1_dict["y_mw"]
    pc2_dict=Pure_component2(num_reactants_condensed,y_mw,R_gas,temp)
    merge!(param_dict,pc1_dict,pc2_dict)
    param_dict["num_reactants_condensed"]=num_reactants_condensed
    println("Generating initial size distribution")
    N_perbin,xs=lognormal(num_bins,total_conc,meansize,size_std,lowersize,uppersize)
    param_dict["N_perbin"]=N_perbin
    
    println("Calculating Dry Core Properties")
    y_core=(4.0/3.0)*pi*((xs*1.0e-6).^3.0) #4/3*pi*radius^3
    y_core=y_core.*core_density_array #mass per particle [kg]
    y_core=y_core./(core_mw*1.0e-3) #moles per particle, changing mw from g/mol to kg/mol
    y_core=y_core*NA #molecules per particle
    y_core=y_core.*N_perbin #molecules/cc representing each size range
    #Calculate a core mass based on the above information [converting from molecules/cc to micrograms/m3]    
    core_mass_array=y_core./NA.*core_mw
    println("Dry core mass = ", sum(core_mass_array)*1E12)
    param_dict["y_core"]=y_core
    param_dict["core_mass_array"]=core_mass_array

    println("Configuring initial condensed phase")
    y_cond=zeros(Float64,num_bins*num_reactants_condensed)
    for step=1:length(xs)
        radius=xs[step]
        water_moles=(y_core[step]*core_dissociation)*(RH/(1.0E0-RH))
        y_cond[step*num_reactants_condensed]=water_moles
    end
    return param_dict,reactants2ind,y_cond
end

function read_configure!(filename::String)
    @printf("Reading Config file %s\n",filename)
    open(filename) do f
        for s in readlines(f)
            if (length(s)>2)
                if s[1]!='#'
                    eval(Meta.parse(s))#eval runs in Module scope while include runs in global scope
                end
            end
        end
    end
end

function run_simulation_aerosol(;use_jacobian::Bool,linsolver::Symbol=:Dense)
    read_configure!("Configure_aerosol.jl")
    param_dict,reactants2ind,y_cond=prepare_aerosol()
    num_reactants,num_reactants_condensed=[param_dict[i] for i in ["num_reactants","num_reactants_condensed"]]
    dy_dt_gas_matrix=zeros(Float64,(num_reactants,num_bins))
    dy_dt=zeros(Float64,num_reactants+num_reactants_condensed*num_bins)
    param_dict["dy_dt_gas_matrix"]=dy_dt_gas_matrix
    param_dict["dydt"]=dy_dt
    param_dict["Current_iter"]=0
    param_dict["Simulation_type"]="aerosol"
    y_init=zeros(Float64,num_reactants+num_reactants_condensed*num_bins)
    for (k,v) in reactants_initial_dict
        y_init[reactants2ind[k]]=v*Cfactor#pbb to molcules/cc
    end
    y_init[num_reactants+1:num_reactants+num_bins*num_reactants_condensed]=y_cond[1:num_bins*num_reactants_condensed]
    println("Solving ODE")
    if use_jacobian
        odefun=ODEFunction(dydt_aerosol!; jac=aerosol_jac!)
        prob = ODEProblem{true}(odefun,y_init,tspan,param_dict)
        param_dict["ShowIterPeriod"]=5
    else
        prob = ODEProblem{true}(dydt_aerosol!,y_init,tspan,param_dict)
        param_dict["ShowIterPeriod"]=500
    end
    sol = solve(prob,CVODE_BDF(linear_solver=linsolver),reltol=1e-4,abstol=1.0e-2,
                tstops=0:batch_step:simulation_time,saveat=batch_step,# save_everystep=true,
                dt=1.0e-6, #Initial step-size
                dtmax=100.0,
                max_order = 5,
                max_convergence_failures = 1000,
                #callback=PositiveDomain(y_init,abstol=1.0e-2)
                #isoutofdomain=(u,p,t) -> any(x -> x < 0, u)
                )
    sol_mtx=transpose(sol)
    aerosol_mtx=sol_mtx[1:end,num_reactants+1:num_reactants+num_bins*num_reactants_condensed]
    t_length=size(aerosol_mtx)[1]
    mw_array=param_dict["y_mw"]
    SOA_array=[sum((sum(reshape(aerosol_mtx[i,1:end],(num_reactants_condensed,num_bins))
                               ,dims=2).*mw_array./NA)[1:end-1]#exclude H2O at the end
                  ) for i in 1:t_length]*1E12

    return sol,reactants2ind,SOA_array,num_reactants,param_dict
end

function run_simulation_aerosol_adjoint(;linsolver::Symbol=:Dense)
    #read_configure!("Configure_aerosol.jl")
    if isfile("/data/aerosol_sol.store")
        read_configure!("Configure_aerosol.jl")
        odefun=ODEFunction(dydt_aerosol!; jac=aerosol_jac!)
        println("Found caching of aerosol simulation")
        sol,param_dict=deserialize("/data/aerosol_sol.store")
    else
        println("No caching, start aerosol simulation")
        sol,_,_,_,param_dict=run_simulation_aerosol(use_jacobian=true,linsolver=linsolver)
        println("Caching solution")
        serialize("/data/aerosol_sol.store",(sol,param_dict))
    end
    num_reactants,num_reactants_condensed,num_eqns=[param_dict[i] for i in ["num_reactants","num_reactants_condensed","num_eqns"]]
    println("Preparing Adjoint Problem")
    t0,tF=tspan
    tspan_adj=(tF,t0)
    len_y=num_reactants+num_bins*num_reactants_condensed
    mw_array=param_dict["y_mw"]
    lambda_init=zeros(Float64,(1,len_y))#DiffEq.jl version seems incorrect
    SOA_mass_jac!(lambda_init,mw_array,NA,num_reactants,num_reactants_condensed,num_bins)#adopting KPP paper I
    #println(lambda_init)
    param_dict["sol"]=sol
    param_dict["jac_mtx"]=zeros(Float64,(len_y,len_y))
    param_dict["Current_iter"]=0
    param_dict["ShowIterPeriod"]=5
    param_dict["Simulation_type"]="adjoint"
    param_dict["jac_cache"]=DiffEqDiffTools.JacobianCache(zeros(Float64,len_y),Val{:forward},Float64,Val{true})
    odefun_adj=ODEFunction(sensitivity_adjoint_dldt!,jac=sensitivity_adjoint_jac!)
    prob_adj=ODEProblem{true}(odefun_adj,lambda_init,tspan_adj,param_dict)
    println("Solving Adjoint Problem")
    lambda_sol=solve(prob_adj,CVODE_BDF(linear_solver=:Dense),reltol=1e-8,abstol=1e-5,
                     tstops=simulation_time:-batch_step:0.,saveat=-batch_step,
                     dt=-1e-8,dtmax=50.0,max_order=5,max_convergence_failures=1000)
    println("Preparing Integration")
    tstops=[t for t in 0:batch_step:simulation_time]
    num_tstops=length(tstops)
    stoich_mtx=param_dict["stoich_mtx"]
    stoich_list=param_dict["stoich_list"]
    reactants_list=param_dict["reactants_list"]
    dSOA_mass_drate=zeros(Float64,(num_eqns,num_tstops))
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
    end
    return dSOA_mass_drate
end

function sensitivity_mtx2dSOA(S,t::Real,integrator)
    p=integrator.p
    mw_array,num_reactants,num_reactants_condensed,num_eqns=[p[i] for i in ["y_mw","num_reactants","num_reactants_condensed","num_eqns"]]
    y_len=num_reactants+num_bins*num_reactants_condensed
    dSOA_dy=zeros(Float64,(1,y_len))
    SOA_mass_jac!(dSOA_dy,mw_array,NA,num_reactants,num_reactants_condensed,num_bins)
    #println(dSOA_dy)
    println(size(S))
    println(S[1:100])
    return reshape(dSOA_dy * reshape(S,(y_len,num_eqns)),num_eqns)
end

function run_simulation_gas()
    read_configure!("Configure_gas.jl")
    param_dict,reactants2ind=prepare_gas()
    num_reactants=param_dict["num_reactants"]
    reactants_initial=zeros(Float64,num_reactants)
    for (k,v) in reactants_initial_dict
        reactants_initial[reactants2ind[k]]=v*Cfactor#pbb to molcules/cc
    end
    println("Solving ODE")
    param_dict["Current_iter"]=0
    param_dict["ShowIterPeriod"]=100
    param_dict["Simulation_type"]="gas"
    #odefun=ODEFunction(dydt!; jac=gas_jac!)
    prob = ODEProblem{true}(dydt!,reactants_initial,tspan,
                            param_dict,
                            #(dy,rate_values,J,stoich_mtx,stoich_list,reactants_list,RO2_inds,num_eqns,num_reactants)
                            )
    sol = solve(prob,CVODE_BDF(linear_solver=:Dense),reltol=1e-6,abstol=1.0e-3,
                tstops=0:batch_step:simulation_time,saveat=batch_step,# save_everystep=true,
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
