function SOA_mass(y,mw_array,num_reactants,num_reactants_condensed,num_bins)
    aerosol_mtx=reshape(y[num_reactants+1:num_reactants+num_reactants_condensed*num_bins],(num_reactants_condensed,num_bins))
    return (sum(aerosol_mtx,dims=2).*mw_array./NA)[1:num_reactants_condensed-1]*1E12
end

function SOA_mass_jac!(dSOA_dy,mw_array,num_reactants,num_reactants_condensed,num_bins)
    for bin_ind in 1:num_bins
        start_ind=num_reactants+1+(bin_ind-1)*num_reactants_condensed
        stop_ind=num_reactants+bin_ind*num_reactants_condensed-1#exclude water in the end
        dSOA_dy[1,start_ind:stop_ind]=((1.0).*mw_array./NA)[1:num_reactants_condensed-1]*1E12
    end
    nothing
end

function loss_gain_drate_values!(num_reactants::Int,num_eqns::Int,
                       reactants::Array{Float64,1},#num_reactants
                       stoich_mtx::SparseMatrixCSC{Float64,Int64},#num_reactants*num_eqns
                       stoich_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, both reac and prod
                       reactants_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, only reac
                       lossgain_jac_mtx#num_output(dydt)*num_input(rate_values) num_reactants*num_eqns 
                       )
    for eqn_ind in 1:num_eqns
        num_reacs,stoichvec,indvec=reactants_list[eqn_ind]
        num_stoichs,_,stoich_indvec=stoich_list[eqn_ind]
        for y_ind in 1:num_reacs
            prod=1.
            for i in 1:num_reacs
                reactant_ind=indvec[i]
                stoich=stoichvec[i]
                prod*=reactants[reactant_ind]^stoich#reactants_list come from reactants_mtx (for catalyse A+B=A+C)
            end
            for i in 1:num_stoichs
                reactant_ind=stoich_indvec[i]
                lossgain_jac_mtx[reactant_ind,eqn_ind]+=stoich_mtx[reactant_ind,eqn_ind]*prod*(-1)
            end
        end
    end
    return lossgain_jac_mtx
end