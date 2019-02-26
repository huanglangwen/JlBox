function Partition!(y::Array{Float64,1},dy_dt::Array{Float64,1},dy_dt_gas_matrix::Array{Float64,2},C_g_i_t::Array{Float64,1},
                    num_bins::Integer,num_reactants::Integer,num_reactants_condensed::Integer,include_inds::Array{Integer,1},
                    mw_array,density_input,gamma_gas,alpha_d_org,DStar_org,Psat,N_perbin::Array{Float64,1},
                    core_diss::Real,y_core::Array{Float64,1},core_mass_array::Array{Float64,1},core_density_array::Array{Float64,1},
                    NA::Real,sigma::Real,R_gas::Real,Model_temp::Real)
    size_array=zeros(Float64,num_bins)
    total_SOA_mass_array=zeros(Float64,num_bins)
    mass_array=zeros(Float64,num_reactants_condensed+1)
    density_array=zeros(Float64,num_reactants_condensed+1)
    #dy_dt_gas_matrix_sum=zeros(Float64,num_reactants)
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
        #println("k_i_m_t,C_g_i_t,Cstar_i_m_t,dm_dt")
        #println(k_i_m_t[end],",",C_g_i_t[end],",",Cstar_i_m_t[end],",",dm_dt[end])
        #println("Pressure_eq,kelvin_factor,y_mole_frac,Psat")
        #println(Pressure_eq[end],",",kelvin_factor[end],",",y_mole_fractions[end],",",Psat[end])

        #ASSIGN dy_dt_gas_matrix
        for ind=1:length(include_inds)
            #dy_dt_gas_matrix_sum[include_inds[ind]]+=dm_dt[ind]
            dy_dt_gas_matrix[include_inds[ind],size_step]=dm_dt[ind]
        end

        dy_dt[start_ind:stop_ind]=dm_dt
    end
    dy_dt[1:num_reactants]=dy_dt[1:num_reactants].-sum(dy_dt_gas_matrix,dims=2)#dy_dt_gas_matrix_sum
    total_SOA_mass=sum(total_SOA_mass_array)*1.0E12
    return dy_dt,total_SOA_mass
end
