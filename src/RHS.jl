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
    if config isa GasConfig
        p["Current_iter"]+=1
        citer=p["Current_iter"]
        if citer%(p["ShowIterPeriod"])==0
            println(config.io, "Current Iteration: $(citer), time_step: $(t)")
            flush(config.io)
        end
    end
    nothing#return dydt
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
        config.sigma,config.temp)
    if config isa AerosolConfig
        p["Current_iter"]+=1
        citer=p["Current_iter"]
        if citer%(p["ShowIterPeriod"])==0
            println(config.io, "Current Iteration: $(citer), time_step: $(t), SOA(ug/m3): $(total_SOA_mass)")
            flush(config.io)
        end
    end
    nothing#return dy_dt
end

function sensitivity_adjoint_dldt!(dldt,lambda,p,t)
    jacobian_from_sol!(p,t)#jacobian_from_sol!(p,t)
    jac_mtx=p["jac_mtx"]
    dldt.= reshape(- lambda' * jac_mtx, :)#adopting KPP paper I
    p["Current_iter"]+=1
    citer=p["Current_iter"]
    adjointconfig = p["adjointconfig"]
    if citer%(p["ShowIterPeriod"])==0
        num_reactants=p["num_reactants"]
        println(adjointconfig.io, "Current Iteration: $(citer), time_step: $(t), sum(lambda_gas): $(sum(lambda[1:num_reactants])), sum(dldt_gas): $(sum(dldt[1:num_reactants])), sum(lambda): $(sum(lambda))")
    end
    nothing
end