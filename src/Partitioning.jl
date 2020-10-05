macro partition_kernel() 
    return esc(quote
        total_moles=sum(temp_array)+y_core[size_step]*core_diss
        y_mole_fractions=temp_array./total_moles

        mass_array[1:num_reactants_condensed]=temp_array.*mw_array./NA
        mass_array[num_reactants_condensed+1]=core_mass_array[size_step]
        density_array[1:num_reactants_condensed]=density_input[1:num_reactants_condensed]
        density_array[num_reactants_condensed+1]=core_density_array[size_step]
        
        #aw_array[size_step]=temp_array[num_reactants_condensed]/total_moles
        total_mass=sum(mass_array)
        total_mass= total_mass > 0 ? total_mass : core_mass_array[size_step]#fix negative size_array^3
        mass_fractions_array=mass_array./total_mass

        #density=1.0/(sum(mass_fractions_array./density_array))
        invdensity=sum(mass_fractions_array./density_array)
        invdensity= invdensity > 0 ? invdensity : 1/core_density_array[size_step]
        
        size_array[size_step]=((3.0*((total_mass*1.0E3)/(N_perbin[size_step]*1.0E6)))*invdensity/(4.0*pi))^(1.0/3.0)

        Kn=gamma_gas./size_array[size_step]
        Inverse_Kn=1.0./Kn
        Correction_part1=(1.33.+0.71*Inverse_Kn)./(1.0.+Inverse_Kn)
        Correction_part2=(4.0*(1.0.-alpha_d_org))./(3.0*alpha_d_org)
        Correction_part3=1.0.+(Correction_part1.+Correction_part2).*Kn
        Correction=1.0./Correction_part3

        kelvin_factor=exp.((4.0*mw_array*1.0E-3*sigma*invdensity)/(R_gas*Model_temp*size_array[size_step]*2.0))
        
        Pressure_eq=kelvin_factor.*y_mole_fractions.*Psat*101325.0

        Cstar_i_m_t=Pressure_eq*(NA/(8.3144598E6*Model_temp))

        k_i_m_t_part1=DStar_org.*Correction
        k_i_m_t=4.0*pi*size_array[size_step]*1.0E2*N_perbin[size_step]*k_i_m_t_part1

        dm_dt=k_i_m_t.*(C_g_i_t-Cstar_i_m_t)
    end)
end

function Partition!(y::Array{FT_DUAL1,1},dy_dt::Array{FT_DUAL2,1},dy_dt_gas_matrix::Array{FT2,2},C_g_i_t::Array{FT_DUAL1,1},
                    num_bins::Integer,num_reactants::Integer,num_reactants_condensed::Integer,include_inds::Array{Integer,1},
                    mw_array,density_input::Array{FT3,1},gamma_gas::Array{FT2,1},alpha_d_org::Array{FT2,1},DStar_org::Array{FT2,1},
                    Psat::Array{FT3,1},N_perbin::Array{FT3,1},core_diss,y_core::Array{FT3,1},core_mass_array::Array{FT3,1},
                    core_density_array::Array{FT3,1},sigma,Model_temp) where {FT_DUAL1<:Real,FT_DUAL2<:Real,FT2<:Real,FT3<:Real}
    size_array=zeros(eltype(y),num_bins)
    total_SOA_mass_array=zeros(eltype(y),num_bins)
    mass_array=zeros(eltype(y),num_reactants_condensed+1)
    density_array=zeros(FT3,num_reactants_condensed+1)
    fill!(dy_dt_gas_matrix,0.)
    #dy_dt_gas_matrix_sum=zeros(Real,num_reactants)
    for size_step=1:num_bins
        start_ind=num_reactants+1+((size_step-1)*num_reactants_condensed)
        stop_ind=num_reactants+(size_step*num_reactants_condensed)
        temp_array=y[start_ind:stop_ind]
        @partition_kernel()
        total_SOA_mass_array[size_step]=sum(mass_array[1:num_reactants_condensed-1])

        #ASSIGN dy_dt_gas_matrix
        @inbounds for ind=1:length(include_inds)
            #dy_dt_gas_matrix_sum[include_inds[ind]]+=dm_dt[ind]
            dy_dt_gas_matrix[include_inds[ind],size_step]=dm_dt[ind]
        end

        dy_dt[start_ind:stop_ind]=dm_dt
    end
    dy_dt[1:num_reactants]=dy_dt[1:num_reactants].-sum(dy_dt_gas_matrix,dims=2)#dy_dt_gas_matrix_sum
    total_SOA_mass=sum(total_SOA_mass_array)*1.0E12
    return dy_dt,total_SOA_mass
end
