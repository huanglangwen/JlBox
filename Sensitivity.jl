function SOA_mass(y,mw_array,NA,num_reactants,num_reactants_condensed,num_bins)
    aerosol_mtx=reshape(y[num_reactants+1:num_reactants+num_reactants_condensed*num_bins],(num_reactants_condensed,num_bins))
    return (sum(aerosol_mtx,dims=2).*mw_array./NA)[1:num_reactants_condensed-1]*1E12
end

function SOA_mass_jac(dSOA_dy,y,mw_array,NA,num_reactants,num_reactants_condensed,num_bins)
    for bin_ind in 1:num_bins
        start_ind=num_reactants+1+(bin_ind-1)*num_reactants_condensed
        stop_ind=num_reactants+bin_ind*num_reactants_condensed-1#exclude water in the end
        dSOA_dy[start_ind:stop_ind]=((1.0).*mw_array./NA)[1:num_reactants_condensed-1]
    end
    nothing
end