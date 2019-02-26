using StaticArrays
using SparseArrays
function loss_gain_jac!(num_reactants::Int,num_eqns::Int,
                       reactants::Array{Float64,1},#num_reactants
                       stoich_mtx::SparseMatrixCSC{Float64,Int64},#num_reactants*num_eqns
                       stoich_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, both reac and prod
                       reactants_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, only reac
                       rate_values::Array{Float64,1},#num_eqns
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
    return lossgain_jac_mtx
end

function gas_jac!(jac_mtx,reactants::Array{Float64,1},p::Dict,t::Real)
    rate_values,stoich_mtx,stoich_list,reactants_list,num_eqns,num_reactants=
        [p[ind] for ind in 
            ["rate_values","stoich_mtx","stoich_list","reactants_list",
             "num_eqns","num_reactants"]
        ]
    #Probably have to re-eval rate_values again
    loss_gain_jac!(num_reactants,num_eqns,reactants,stoich_mtx,stoich_list,reactants_list,rate_values,jac_mtx)
end

function Partition_jac!(y_jac,y::Array{Float64,1},dy_dt::Array{Float64,1},dy_dt_gas_matrix::Array{Float64,2},C_g_i_t::Array{Float64,1},
                        num_bins::Integer,num_reactants::Integer,num_reactants_condensed::Integer,include_inds::Array{Integer,1},
                        mw_array,density_input,gamma_gas,alpha_d_org,DStar_org,Psat,N_perbin::Array{Float64,1},
                        core_diss::Real,y_core::Array{Float64,1},core_mass_array::Array{Float64,1},core_density_array::Array{Float64,1},
                        NA::Real,sigma::Real,R_gas::Real,Model_temp::Real)
    #y_jac: Jacobian matrix (num_output*num_input) num_output==num_input==num_reactants+num_bins*num_reactants_condensed
    size_array=zeros(Float64,num_bins)
    total_SOA_mass_array=zeros(Float64,num_bins)
    mass_array=zeros(Float64,num_reactants_condensed+1)
    density_array=zeros(Float64,num_reactants_condensed+1)
    for size_step=1:num_bins
        start_ind=num_reactants+1+((size_step-1)*num_reactants_condensed)
        stop_ind=num_reactants+(size_step*num_reactants_condensed)
        temp_array=y[start_ind:stop_ind]
        total_moles=sum(temp_array)+y_core[size_step]*core_diss
        y_mole_fractions=temp_array./total_moles

        mass_array[1:num_reactants_condensed]=temp_array.*mw_array./NA
        mass_array[num_reactants_condensed+1]=core_mass_array[size_step]
        density_array[1:num_reactants_condensed]=density_input[1:num_reactants_condensed]
        density_array[num_reactants_condensed+1]=core_density_array[size_step]
        
        total_SOA_mass_array[size_step]=sum(mass_array[1:num_reactants_condensed-1])
        #aw_array[size_step]=temp_array[num_reactants_condensed]/total_moles
        total_mass=sum(mass_array)
        mass_fractions_array=mass_array./total_mass

        density=1.0/(sum(mass_fractions_array./density_array))
        
        size_array[size_step]=((3.0*((total_mass*1.0E3)/(N_perbin[size_step]*1.0E6)))/(4.0*pi*density))^(1.0/3.0)

        Kn=gamma_gas./size_array[size_step]
        Inverse_Kn=1.0./Kn
        Correction_part1=(1.33.+0.71*Inverse_Kn)./(1.0.+Inverse_Kn)
        Correction_part2=(4.0*(1.0.-alpha_d_org))./(3.0*alpha_d_org)
        Correction_part3=1.0.+(Correction_part1.+Correction_part2).*Kn
        Correction=1.0./Correction_part3

        kelvin_factor=exp.((4.0*mw_array*1.0E-3*sigma)/(R_gas*Model_temp*size_array[size_step]*2.0*density))
        
        Pressure_eq=kelvin_factor.*y_mole_fractions.*Psat*101325.0

        Cstar_i_m_t=Pressure_eq*(NA/(8.3144598E6*Model_temp))

        k_i_m_t_part1=DStar_org.*Correction
        k_i_m_t=4.0*pi*size_array[size_step]*1.0E2*N_perbin[size_step]*k_i_m_t_part1

        dm_dt=k_i_m_t.*(C_g_i_t-Cstar_i_m_t)

        #ASSIGN dy_dt_gas_matrix
        for ind=1:length(include_inds)
            dy_dt_gas_matrix[include_inds[ind],size_step]=dm_dt[ind]
        end

        dy_dt[start_ind:stop_ind]=dm_dt
    end
    dy_dt[1:num_reactants]=dy_dt[1:num_reactants]-sum(dy_dt_gas_matrix,dims=2)
    total_SOA_mass=sum(total_SOA_mass_array)*1.0E12
    return dy_dt,total_SOA_mass
end