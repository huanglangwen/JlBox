function select_jacobian(diff_method, len_y)
    if diff_method == "finite"
        jac_cache = FiniteDiff.JacobianCache(zeros(Float64,len_y),zeros(Float64,len_y),Val{:forward},Float64)
        jac! = (jac_mtx,y,p,t,fillzero::Bool=true)-> begin
            if fillzero
                fill!(jac_mtx,0.0)
            end
            FiniteDiff.finite_difference_jacobian!(jac_mtx,(dydt,y)->dydt_aerosol!(dydt,y,p,t),y,jac_cache)
        end
    elseif diff_method == "coarse_seeding"
        jac! = aerosol_jac_coarse_seeding!
    elseif diff_method == "fine_seeding"
        jac! = aerosol_jac_fine_seeding!
    elseif diff_method == "fine_analytical"
        jac! = aerosol_jac_fine_analytical!
    elseif diff_method == "coarse_analytical"
        jac! = aerosol_jac_coarse_analytical!
    elseif diff_method == "gas"
        jac! = gas_jac!
    else
        println("ERROR: can't recognize diff type: ",diff_method)
        exit(-1)
    end
    jac!
end

function loss_gain_jac!(num_reactants::Int,num_eqns::Int,
                       reactants::Array{<:Real,1},#num_reactants
                       stoich_mtx::SparseMatrixCSC{Float64,Int64},#num_reactants*num_eqns
                       stoich_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, both reac and prod
                       reactants_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, only reac
                       rate_values::Array{<:Real,1},#num_eqns
                       lossgain_jac_mtx#::Array{Float64,2}#SparseMatrixCSC{Float64,Int64},#num_output(dydt)*num_input(y)
                       )
    #lossgain_jac_mtx=spzeros(num_reactants,num_reactants)#num_output(dydt)*num_input(y)
    for eqn_ind in 1:num_eqns
        num_reacs,stoichvec,indvec=reactants_list[eqn_ind]
        num_stoichs,_,stoich_indvec=stoich_list[eqn_ind]
        for y_ind in 1:num_reacs
            prod=rate_values[eqn_ind]
            for i in [i for i in 1:num_reacs if i!=y_ind]
                reactant_ind=indvec[i]
                stoich=stoichvec[i]
                prod*=reactants[reactant_ind]^stoich#reactants_list come from reactants_mtx (for catalyse A+B=A+C)
            end
            reactant_y_ind=indvec[y_ind]
            stoich_y::Integer=stoichvec[y_ind]
            prod*=stoich_y*reactants[reactant_y_ind]^(stoich_y-1)
            for i in 1:num_stoichs
                reactant_ind=stoich_indvec[i]
                lossgain_jac_mtx[reactant_ind,reactant_y_ind]+=stoich_mtx[reactant_ind,eqn_ind]*prod*(-1)
            end
        end
    end
    nothing
end

function gas_jac!(jac_mtx,reactants::Array{<:Real,1},p::Dict,t::Real,fillzero::Bool=true)
    rate_values,stoich_mtx,stoich_list,reactants_list,num_eqns,num_reactants=
        [p[ind] for ind in 
            ["rate_values","stoich_mtx","stoich_list","reactants_list",
             "num_eqns","num_reactants"]
        ]
    #Probably have to re-eval rate_values again
    if fillzero
        fill!(jac_mtx,0.0)
    end
    loss_gain_jac!(num_reactants,num_eqns,reactants,stoich_mtx,stoich_list,reactants_list,rate_values,jac_mtx)
end

function loss_gain_jac_v!(num_reactants::Int,num_eqns::Int,
                       reactants::Array{<:Real,1},#num_reactants
                       stoich_mtx::SparseMatrixCSC{Float64,Int64},#num_reactants*num_eqns
                       stoich_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, both reac and prod
                       reactants_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, only reac
                       rate_values::Array{<:Real,1},#num_eqns
                       v::Array{Float64,1},#num_reactants
                       lossgain_jac_v::Array{Float64,1}#SparseMatrixCSC{Float64,Int64},#num_output(dydt)
                       )
    #lossgain_jac_mtx=spzeros(num_reactants,num_reactants)#num_output(dydt)*num_input(y)
    for eqn_ind in 1:num_eqns
        num_reacs,stoichvec,indvec=reactants_list[eqn_ind]
        num_stoichs,_,stoich_indvec=stoich_list[eqn_ind]
        for y_ind in 1:num_reacs
            prod=rate_values[eqn_ind]
            for i in [i for i in 1:num_reacs if i!=y_ind]
                reactant_ind=indvec[i]
                stoich=stoichvec[i]
                prod*=reactants[reactant_ind]^stoich#reactants_list come from reactants_mtx (for catalyse A+B=A+C)
            end
            reactant_y_ind=indvec[y_ind]
            stoich_y::Integer=stoichvec[y_ind]
            prod*=stoich_y*reactants[reactant_y_ind]^(stoich_y-1)
            for i in 1:num_stoichs
                reactant_ind=stoich_indvec[i]
                lossgain_jac_v[reactant_ind]+=v[reactant_y_ind]*stoich_mtx[reactant_ind,eqn_ind]*prod*(-1)
            end
        end
    end
    nothing
end

function gas_jac_v!(jac_v::Array{Float64,1},v::Array{Float64,1},reactants::Array{Float64,1},p::Dict,t::Float64,fillzero::Bool=true)
    rate_values,stoich_mtx,stoich_list,reactants_list,num_eqns,num_reactants=
        [p[ind] for ind in 
            ["rate_values","stoich_mtx","stoich_list","reactants_list",
             "num_eqns","num_reactants"]
        ]
    #Probably have to re-eval rate_values again
    if fillzero
        fill!(jac_v,0.0)
    end
    loss_gain_jac_v!(num_reactants,num_eqns,reactants,stoich_mtx,stoich_list,reactants_list,rate_values,v,jac_v)
end

function Partition_jac!(y_jac,y::Array{Float64,1},C_g_i_t::Array{Float64,1},
                        num_bins::Integer,num_reactants::Integer,num_reactants_condensed::Integer,include_inds::Array{Integer,1},
                        mw_array,density_input,gamma_gas,alpha_d_org,DStar_org,Psat,N_perbin::Array{Float64,1},
                        core_diss::Real,y_core::Array{Float64,1},core_mass_array::Array{Float64,1},core_density_array::Array{Float64,1},
                        sigma::Real,Model_temp::Real)
    #y_jac: Jacobian matrix (num_output*num_input) num_output==num_input==num_reactants+num_bins*num_reactants_condensed
    size_array=zeros(Float64,num_bins)
    #total_SOA_mass_array=zeros(Float64,num_bins)
    mass_array=zeros(Float64,num_reactants_condensed+1)
    density_array=zeros(Float64,num_reactants_condensed+1)
    #DC_g_i_t=spzeros(num_reactants_condensed,num_reactants)
    #Ddm_dt_Dy_gas_sum=spzeros(num_reactants_condensed,num_reactants)
    DC_g_i_t=zeros(num_reactants_condensed,num_reactants)
    Ddm_dt_Dy_gas_sum=zeros(num_reactants_condensed,num_reactants)
    for i in 1:num_reactants_condensed
        DC_g_i_t[i,include_inds[i]]=1
    end
    for size_step=1:num_bins
        start_ind=num_reactants+1+((size_step-1)*num_reactants_condensed)
        stop_ind=num_reactants+(size_step*num_reactants_condensed)
        temp_array=y[start_ind:stop_ind]
        @partition_kernel()

        #ASSIGN dy_dt_gas_matrix
        #for ind=1:length(include_inds)
        #    dy_dt_gas_matrix[include_inds[ind],size_step]=dm_dt[ind]
        #end

        #dy_dt[start_ind:stop_ind]=dm_dt

        #=================Jacobian, Input=y_gas======================#
        begin
            Ddm_dt_Dy_gas=k_i_m_t.*DC_g_i_t
            y_jac[start_ind:stop_ind,1:num_reactants]=Ddm_dt_Dy_gas#num_condensed*num_reactants
            Ddm_dt_Dy_gas_sum.+=Ddm_dt_Dy_gas
        end
        #=================Jacobian, Input=y_bins======================#
        begin
            #Dtemp_array=sparse(1:num_reactants_condensed,1:num_reactants_condensed,ones(num_reactants_condensed))#I, num_condensed*num_condensed
            Dtemp_array=Matrix{Float64}(I,num_reactants_condensed,num_reactants_condensed)
            Dtotal_moles=ones(1,num_reactants_condensed)
            Dy_mole_fractions=1/total_moles.*Dtemp_array.-(temp_array./total_moles^2).*Dtotal_moles#num_condensed*num_condensed
            #Dmass_array=spzeros(num_reactants_condensed+1,num_reactants_condensed)
            Dmass_array=zeros(num_reactants_condensed+1,num_reactants_condensed)
            Dmass_array[1:num_reactants_condensed,:]=(mw_array./NA).*Dtemp_array
            Dtotal_mass=sum(Dmass_array,dims=1)#1*num_condensed
            Dmass_fractions_array=Dmass_array./total_mass-(mass_array./(total_mass^2)).*Dtotal_mass#(num_condensed+1)*num_condensed
            density=1/invdensity
            Ddensity=-(density^2).*sum(Dmass_fractions_array./density_array,dims=1)#1*num_condensed
            Dsize=1/3*size_array[size_step]^(-2)*3.0*1E3/(N_perbin[size_step]*1E6*4*pi)*(Dtotal_mass./density.-total_mass/(density^2).*Ddensity)#1*num_condensed
            DKn=-gamma_gas./(size_array[size_step]^2).*Dsize#num_condensed*num_condensed
            DCorrection_part1=(1.33-0.71)./((Kn.+1).^2).*DKn#num_condensed*num_condensed
            DCorrection_part3=Kn.*DCorrection_part1+Correction_part1.*DKn#num_condensed*num_condensed
            DCorrection=-(Correction.^2).*DCorrection_part3#num_condensed*num_condensed
            Dkelvin_factor=(kelvin_factor.*(4.0*mw_array*1.0E-3*sigma)./(R_gas*Model_temp*2.0)*(-1)./(size_array[size_step]*density)^2).*(density*Dsize+size_array[size_step]*Ddensity)#num_condensed*num_condensed
            DPressure_eq=(Psat*101325.0).*(y_mole_fractions.*Dkelvin_factor+kelvin_factor.*Dy_mole_fractions)#num_condensed*num_condensed
            DCstar_i_m_t=(NA/(8.3144598E6*Model_temp)).*DPressure_eq#num_condensed*num_condensed
            Dk_i_m_t_part1=DStar_org.*DCorrection#num_condensed*num_condensed
            Dk_i_m_t=4.0*pi*1.0E2*N_perbin[size_step].*(size_array[size_step].*Dk_i_m_t_part1.+k_i_m_t_part1.*Dsize)#num_condensed*num_condensed
            Ddm_dt=(C_g_i_t-Cstar_i_m_t).*Dk_i_m_t-k_i_m_t.*DCstar_i_m_t#num_condensed*num_condensed
            y_jac[start_ind:stop_ind,start_ind:stop_ind]=Ddm_dt#num_condensed*num_condensed
            y_jac[1:num_reactants,start_ind:stop_ind]=-transpose(DC_g_i_t)*Ddm_dt#num_reactants*num_condensed
        end
    end
    #dy_dt[1:num_reactants]=dy_dt[1:num_reactants]-sum(dy_dt_gas_matrix,dims=2)
    #total_SOA_mass=sum(total_SOA_mass_array)*1.0E12

    #Jacobian
    Ddy_dt_gas_matrix_sum_Dy_gas=transpose(DC_g_i_t)*Ddm_dt_Dy_gas_sum#num_reactants*num_reactants
    y_jac[1:num_reactants,1:num_reactants].-=Ddy_dt_gas_matrix_sum_Dy_gas

    nothing
    #return dy_dt,total_SOA_mass
end

function Partition_jac_simplified!(y_jac,y::Array{Float64,1},C_g_i_t::Array{Float64,1},
                        num_bins::Integer,num_reactants::Integer,num_reactants_condensed::Integer,include_inds::Array{Integer,1},
                        mw_array,density_input,gamma_gas,alpha_d_org,DStar_org,Psat,N_perbin::Array{Float64,1},
                        core_diss::Real,y_core::Array{Float64,1},core_mass_array::Array{Float64,1},core_density_array::Array{Float64,1},
                        sigma::Real,Model_temp::Real)
    #y_jac: Jacobian matrix (num_output*num_input) num_output==num_input==num_reactants+num_bins*num_reactants_condensed
    size_array=zeros(Float64,num_bins)
    #total_SOA_mass_array=zeros(Float64,num_bins)
    mass_array=zeros(Float64,num_reactants_condensed+1)
    density_array=zeros(Float64,num_reactants_condensed+1)
    #DC_g_i_t=spzeros(num_reactants_condensed,num_reactants)
    #Ddm_dt_Dy_gas_sum=spzeros(num_reactants_condensed,num_reactants)
    DC_g_i_t=zeros(num_reactants_condensed,num_reactants)
    Ddm_dt_Dy_gas_sum=zeros(num_reactants_condensed,num_reactants)
    for i in 1:num_reactants_condensed
        DC_g_i_t[i,include_inds[i]]=1
    end
    for size_step=1:num_bins
        start_ind=num_reactants+1+((size_step-1)*num_reactants_condensed)
        stop_ind=num_reactants+(size_step*num_reactants_condensed)
        temp_array=y[start_ind:stop_ind]
        @partition_kernel()

        #ASSIGN dy_dt_gas_matrix
        #for ind=1:length(include_inds)
        #    dy_dt_gas_matrix[include_inds[ind],size_step]=dm_dt[ind]
        #end

        #dy_dt[start_ind:stop_ind]=dm_dt

        #=================Jacobian, Input=y_gas======================#
        begin
            Ddm_dt_Dy_gas=k_i_m_t.*DC_g_i_t
            y_jac[start_ind:stop_ind,1:num_reactants]=Ddm_dt_Dy_gas#num_condensed*num_reactants
            Ddm_dt_Dy_gas_sum.+=Ddm_dt_Dy_gas
        end
        #=================Jacobian, Input=y_bins======================#
        begin
            #Dtemp_array=sparse(1:num_reactants_condensed,1:num_reactants_condensed,ones(num_reactants_condensed))#I, num_condensed*num_condensed
            Dtemp_array=Matrix{Float64}(I,num_reactants_condensed,num_reactants_condensed)
            Dtotal_moles=ones(1,num_reactants_condensed)
            Dy_mole_fractions=1/total_moles.*Dtemp_array.-(temp_array./total_moles^2).*Dtotal_moles#num_condensed*num_condensed
            #Dmass_array=spzeros(num_reactants_condensed+1,num_reactants_condensed)
            #Dmass_array=zeros(num_reactants_condensed+1,num_reactants_condensed)
            #Dmass_array[1:num_reactants_condensed,:]=(mw_array./NA).*Dtemp_array
            #Dtotal_mass=sum(Dmass_array,dims=1)#1*num_condensed
            #Dmass_fractions_array=Dmass_array./total_mass-(mass_array./(total_mass^2)).*Dtotal_mass#(num_condensed+1)*num_condensed
            #Ddensity=-(density^2).*sum(Dmass_fractions_array./density_array,dims=1)#1*num_condensed
            #Dsize=1/3*size_array[size_step]^(-2)*3.0*1E3/(N_perbin[size_step]*1E6*4*pi)*(Dtotal_mass./density.-total_mass/(density^2).*Ddensity)#1*num_condensed
            #DKn=-gamma_gas./(size_array[size_step]^2).*Dsize#num_condensed*num_condensed
            #DCorrection_part1=(1.33-0.71)./((Kn.+1).^2).*DKn#num_condensed*num_condensed
            #DCorrection_part3=Kn.*DCorrection_part1+Correction_part1.*DKn#num_condensed*num_condensed
            #DCorrection=-(Correction.^2).*DCorrection_part3#num_condensed*num_condensed
            #Dkelvin_factor=(kelvin_factor.*(4.0*mw_array*1.0E-3*sigma)./(R_gas*Model_temp*2.0)*(-1)./(size_array[size_step]*density)^2).*(density*Dsize+size_array[size_step]*Ddensity)#num_condensed*num_condensed
            #DPressure_eq=(Psat*101325.0).*(y_mole_fractions.*Dkelvin_factor+kelvin_factor.*Dy_mole_fractions)#num_condensed*num_condensed
            DPressure_eq=(Psat*101325.0).*kelvin_factor.*Dy_mole_fractions#num_condensed*num_condensed
            DCstar_i_m_t=(NA/(8.3144598E6*Model_temp)).*DPressure_eq#num_condensed*num_condensed
            #Dk_i_m_t_part1=DStar_org.*DCorrection#num_condensed*num_condensed
            #Dk_i_m_t=4.0*pi*1.0E2*N_perbin[size_step].*(size_array[size_step].*Dk_i_m_t_part1.+k_i_m_t_part1.*Dsize)#num_condensed*num_condensed
            #Ddm_dt=(C_g_i_t-Cstar_i_m_t).*Dk_i_m_t-k_i_m_t.*DCstar_i_m_t#num_condensed*num_condensed
            Ddm_dt=-k_i_m_t.*DCstar_i_m_t#num_condensed*num_condensed
            y_jac[start_ind:stop_ind,start_ind:stop_ind].=Ddm_dt#num_condensed*num_condensed
            y_jac[1:num_reactants,start_ind:stop_ind].=-transpose(DC_g_i_t)*Ddm_dt#num_reactants*num_condensed
        end
    end
    #dy_dt[1:num_reactants]=dy_dt[1:num_reactants]-sum(dy_dt_gas_matrix,dims=2)
    #total_SOA_mass=sum(total_SOA_mass_array)*1.0E12

    #Jacobian
    Ddy_dt_gas_matrix_sum_Dy_gas=transpose(DC_g_i_t)*Ddm_dt_Dy_gas_sum#num_reactants*num_reactants
    y_jac[1:num_reactants,1:num_reactants].-=Ddy_dt_gas_matrix_sum_Dy_gas

    nothing
    #return dy_dt,total_SOA_mass
end

function Partition_jac_AD!(y_jac,y::Array{Float64,1},C_g_i_t::Array{Float64,1},
                        num_bins::Integer,num_reactants::Integer,num_reactants_condensed::Integer,include_inds::Array{Integer,1},
                        mw_array,density_input,gamma_gas,alpha_d_org,DStar_org,Psat,N_perbin::Array{Float64,1},
                        core_diss::Real,y_core::Array{Float64,1},core_mass_array::Array{Float64,1},core_density_array::Array{Float64,1},
                        sigma::Real,Model_temp::Real)
    #y_jac: Jacobian matrix (num_output*num_input) num_output==num_input==num_reactants+num_bins*num_reactants_condensed
    size_array=zeros(Real,num_bins)
    #total_SOA_mass_array=zeros(Float64,num_bins)
    mass_array=zeros(Real,num_reactants_condensed+1)
    density_array=zeros(Real,num_reactants_condensed+1)
    #DC_g_i_t=spzeros(num_reactants_condensed,num_reactants)
    #Ddm_dt_Dy_gas_sum=spzeros(num_reactants_condensed,num_reactants)
    DC_g_i_t=zeros(num_reactants_condensed,num_reactants)
    Ddm_dt_Dy_gas_sum=zeros(num_reactants_condensed,num_reactants)
    Ddm_dt_Dy_bin=zeros(num_reactants_condensed,num_reactants_condensed)
    dm_dt_arr=zeros(num_reactants_condensed)
    for i in 1:num_reactants_condensed
        DC_g_i_t[i,include_inds[i]]=1
    end
    for size_step=1:num_bins
        start_ind=num_reactants+1+((size_step-1)*num_reactants_condensed)
        stop_ind=num_reactants+(size_step*num_reactants_condensed)
        temp_array=y[start_ind:stop_ind]
        #=================Jacobian, Input=y_gas======================#
        begin
            @partition_kernel()
            Ddm_dt_Dy_gas=k_i_m_t.*DC_g_i_t
            y_jac[start_ind:stop_ind,1:num_reactants]=Ddm_dt_Dy_gas#num_condensed*num_reactants
            Ddm_dt_Dy_gas_sum.+=Ddm_dt_Dy_gas
        end
        #=================Jacobian, Input=y_bin======================#
        begin
            function dm_dt_func(dm_dt_arr::Array{<:Real,1}, temp_array::Array{<:Real,1})
                @partition_kernel()#captures temp_array
                dm_dt_arr .= dm_dt
                nothing
            end
            ForwardDiff.jacobian!(Ddm_dt_Dy_bin, dm_dt_func, dm_dt_arr, temp_array)
            y_jac[start_ind:stop_ind,start_ind:stop_ind]=Ddm_dt_Dy_bin#num_condensed*num_condensed
            y_jac[1:num_reactants,start_ind:stop_ind]=-transpose(DC_g_i_t)*Ddm_dt_Dy_bin#num_reactants*num_condensed
        end
    end
    #dy_dt[1:num_reactants]=dy_dt[1:num_reactants]-sum(dy_dt_gas_matrix,dims=2)
    #total_SOA_mass=sum(total_SOA_mass_array)*1.0E12

    #Jacobian
    Ddy_dt_gas_matrix_sum_Dy_gas=transpose(DC_g_i_t)*Ddm_dt_Dy_gas_sum#num_reactants*num_reactants
    y_jac[1:num_reactants,1:num_reactants].-=Ddy_dt_gas_matrix_sum_Dy_gas

    nothing
    #return dy_dt,total_SOA_mass
end

#https://github.com/SciML/DiffEqOperators.jl/blob/master/src/jacvec_operators.jl
struct JacVecTag end
function auto_jacvec!(du, f, x, v,
                 cache1 = ForwardDiff.Dual{JacVecTag}.(x, x), # this won't alias
                 cache2 = similar(cache1))
    cache1 .= ForwardDiff.Dual{JacVecTag}.(x, reshape(v, size(x)))
    f(cache2,cache1)
    du .= vec(ForwardDiff.partials.(cache2, 1))
end

function Partition_jac_v_AD!(y_jac_v,v::Array{Float64,1},y::Array{Float64,1},C_g_i_t::Array{Float64,1},
                        num_bins::Integer,num_reactants::Integer,num_reactants_condensed::Integer,include_inds::Array{Integer,1},
                        mw_array,density_input,gamma_gas,alpha_d_org,DStar_org,Psat,N_perbin::Array{Float64,1},
                        core_diss::Real,y_core::Array{Float64,1},core_mass_array::Array{Float64,1},core_density_array::Array{Float64,1},
                        sigma::Real,Model_temp::Real)
    size_array=zeros(Real,num_bins)
    mass_array=zeros(Real,num_reactants_condensed+1)
    density_array=zeros(Real,num_reactants_condensed+1)
    DC_g_i_t=zeros(num_reactants_condensed,num_reactants)
    Ddm_dt_v_Dy_gas_sum=zeros(num_reactants_condensed)
    Ddm_dt_v_Dy_bin=zeros(num_reactants_condensed)
    dm_dt_arr=zeros(num_reactants_condensed)
    cache1=ForwardDiff.Dual{JacVecTag}.(Ddm_dt_v_Dy_bin, Ddm_dt_v_Dy_bin)
    cache2=similar(cache1)
    for i in 1:num_reactants_condensed
        DC_g_i_t[i,include_inds[i]]=1
    end
    for size_step=1:num_bins
        start_ind=num_reactants+1+((size_step-1)*num_reactants_condensed)
        stop_ind=num_reactants+(size_step*num_reactants_condensed)
        temp_array=y[start_ind:stop_ind]
        #=================Jacobian, Input=y_gas======================#
        begin
            @partition_kernel()
            #k_i_m_t: num_reactants_condensed
            Ddm_dt_v_Dy_gas = k_i_m_t.*(DC_g_i_t*v[1:num_reactants])#num_reactants_condensed
            y_jac_v[start_ind:stop_ind] .+= Ddm_dt_v_Dy_gas
            Ddm_dt_v_Dy_gas_sum .+= Ddm_dt_v_Dy_gas
        end
        #=================Jacobian, Input=y_bin======================#
        begin
            function dm_dt_func(dm_dt_arr::Array{<:Real,1}, temp_array::Array{<:Real,1})
                @partition_kernel()#captures temp_array
                dm_dt_arr .= dm_dt
                nothing
            end
            auto_jacvec!(Ddm_dt_v_Dy_bin,dm_dt_func,temp_array,v[start_ind:stop_ind],cache1,cache2)
            y_jac_v[start_ind:stop_ind] .+= Ddm_dt_v_Dy_bin
            y_jac_v[1:num_reactants] .+= -transpose(DC_g_i_t)*Ddm_dt_v_Dy_bin
        end
    end
    #Jacobian
    y_jac_v[1:num_reactants] .+= -transpose(DC_g_i_t)*Ddm_dt_v_Dy_gas_sum

    nothing
    #return dy_dt,total_SOA_mass
end

function aerosol_jac_fine_analytical!(jac_mtx,y::Array{Float64,1},p::Dict,t::Real,fillzero::Bool=true)
    num_reactants,num_reactants_condensed=[p[i] for i in ["num_reactants","num_reactants_condensed"]]
    rate_values,J,RO2_inds=[p[i] for i in ["rate_values","J","RO2_inds"]]
    evaluate_rates_fun=p["evaluate_rates!"]
    config=p["config"]
    time_of_day_seconds=config.start_time+t
    RO2=sum(y[RO2_inds])
    Base.invokelatest(evaluate_rates_fun,time_of_day_seconds,RO2,config.H2O,config.temp,rate_values,J)
    if fillzero
        fill!(jac_mtx,0.0)
    end
    gas_jac!(jac_mtx,y,p,t,false)
    include_inds,N_perbin=[p[i] for i in ["include_inds","N_perbin"]]
    mw_array,density_array,gamma_gas,alpha_d_org,DStar_org,Psat=[p[i] for i in ["y_mw","y_density_array","gamma_gas","alpha_d_org","DStar_org","Psat"]]
    y_core,core_mass_array=[p[i] for i in ["y_core","core_mass_array"]]
    C_g_i_t=y[include_inds]
    Partition_jac!(jac_mtx,y,C_g_i_t,
        config.num_bins,num_reactants,num_reactants_condensed,include_inds,
        mw_array,density_array,gamma_gas,alpha_d_org,DStar_org,Psat,N_perbin,
        config.core_dissociation,y_core,core_mass_array,config.core_density_array,
        config.sigma,config.temp)
    nothing
end

function aerosol_jac_coarse_analytical!(jac_mtx,y::Array{Float64,1},p::Dict,t::Real,fillzero::Bool=true)
    num_reactants,num_reactants_condensed=[p[i] for i in ["num_reactants","num_reactants_condensed"]]
    rate_values,J,RO2_inds=[p[i] for i in ["rate_values","J","RO2_inds"]]
    evaluate_rates_fun=p["evaluate_rates!"]
    config=p["config"]
    time_of_day_seconds=config.start_time+t
    RO2=sum(y[RO2_inds])
    Base.invokelatest(evaluate_rates_fun,time_of_day_seconds,RO2,config.H2O,config.temp,rate_values,J)
    if fillzero
        fill!(jac_mtx,0.0)
    end
    gas_jac!(jac_mtx,y,p,t,false)
    include_inds,N_perbin=[p[i] for i in ["include_inds","N_perbin"]]
    mw_array,density_array,gamma_gas,alpha_d_org,DStar_org,Psat=[p[i] for i in ["y_mw","y_density_array","gamma_gas","alpha_d_org","DStar_org","Psat"]]
    y_core,core_mass_array=[p[i] for i in ["y_core","core_mass_array"]]
    C_g_i_t=y[include_inds]
    Partition_jac_simplified!(jac_mtx,y,C_g_i_t,
        config.num_bins,num_reactants,num_reactants_condensed,include_inds,
        mw_array,density_array,gamma_gas,alpha_d_org,DStar_org,Psat,N_perbin,
        config.core_dissociation,y_core,core_mass_array,config.core_density_array,
        config.sigma,config.temp)
    nothing
end

function aerosol_jac_fine_seeding!(jac_mtx,y::Array{Float64,1},p::Dict,t::Real,fillzero::Bool=true)
    num_reactants,num_reactants_condensed=[p[i] for i in ["num_reactants","num_reactants_condensed"]]
    rate_values,J,RO2_inds=[p[i] for i in ["rate_values","J","RO2_inds"]]
    evaluate_rates_fun=p["evaluate_rates!"]
    config=p["config"]
    time_of_day_seconds=config.start_time+t
    RO2=sum(y[RO2_inds])
    Base.invokelatest(evaluate_rates_fun,time_of_day_seconds,RO2,config.H2O,config.temp,rate_values,J)
    if fillzero
        fill!(jac_mtx,0.0)
    end
    gas_jac!(jac_mtx,y,p,t,false)
    include_inds,N_perbin=[p[i] for i in ["include_inds","N_perbin"]]
    mw_array,density_array,gamma_gas,alpha_d_org,DStar_org,Psat=[p[i] for i in ["y_mw","y_density_array","gamma_gas","alpha_d_org","DStar_org","Psat"]]
    y_core,core_mass_array=[p[i] for i in ["y_core","core_mass_array"]]
    C_g_i_t=y[include_inds]
    Partition_jac_AD!(jac_mtx,y,C_g_i_t,
        config.num_bins,num_reactants,num_reactants_condensed,include_inds,
        mw_array,density_array,gamma_gas,alpha_d_org,DStar_org,Psat,N_perbin,
        config.core_dissociation,y_core,core_mass_array,config.core_density_array,
        config.sigma,config.temp)
    nothing
end

function aerosol_jac_coarse_seeding!(jac_mtx,y::Array{Float64,1},p::Dict,t::Real,fillzero::Bool=true)
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
        config.sigma,config.temp)
    end
    dy_dt=zeros(Real,length(y))
    if fillzero
        fill!(jac_mtx,0.0)
    end
    ForwardDiff.jacobian!(jac_mtx,partition_dydt_fun, dy_dt, y)
    gas_jac!(jac_mtx,y,p,t,false)
    nothing
end

function aerosol_jac_v_fine_seeding!(jac_v::Array{Float64,1},v::Array{Float64,1},y::Array{Float64,1},p::Dict,t::Real,fillzero::Bool=true)
    num_reactants,num_reactants_condensed=[p[i] for i in ["num_reactants","num_reactants_condensed"]]
    rate_values,J,RO2_inds=[p[i] for i in ["rate_values","J","RO2_inds"]]
    evaluate_rates_fun=p["evaluate_rates!"]
    config=p["config"]
    time_of_day_seconds=config.start_time+t
    RO2=sum(y[RO2_inds])
    Base.invokelatest(evaluate_rates_fun,time_of_day_seconds,RO2,config.H2O,config.temp,rate_values,J)
    if fillzero
        fill!(jac_v,0.0)
    end
    gas_jac_v!(jac_v,v,y,p,t,false)
    include_inds,N_perbin=[p[i] for i in ["include_inds","N_perbin"]]
    mw_array,density_array,gamma_gas,alpha_d_org,DStar_org,Psat=[p[i] for i in ["y_mw","y_density_array","gamma_gas","alpha_d_org","DStar_org","Psat"]]
    y_core,core_mass_array=[p[i] for i in ["y_core","core_mass_array"]]
    C_g_i_t=y[include_inds]
    Partition_jac_v_AD!(jac_v,v,y,C_g_i_t,
        config.num_bins,num_reactants,num_reactants_condensed,include_inds,
        mw_array,density_array,gamma_gas,alpha_d_org,DStar_org,Psat,N_perbin,
        config.core_dissociation,y_core,core_mass_array,config.core_density_array,
        config.sigma,config.temp)
    nothing
end

function sensitivity_adjoint_jac!(jac_mtx,lambda,p,t)
    jacobian_from_sol!(p,t)#jacobian_from_sol!(p,t)
    jac_mtx.=(-1).*transpose(p["jac_mtx"])#IMPORTANT jacobian should be the transpose of the original one 
    # since dldt=g(t)-l*J, for ith element in l and jth element in dldt appears at ith line and jth col in the Jacobian matrix
    nothing
end

function jacobian_from_sol!(p::Dict,t::Real)
    sol = p["sol"]
    y = sol(t)
    jac_mtx = p["jac_mtx"]
    jac! = p["jacobian!"]
    jac!(jac_mtx,y,p,t,true)
    nothing
end