
using StaticArrays
using SparseArrays
function mk_reactants_list(num_reactants::Int,num_eqns::Int,
                           reactants_mtx::SparseMatrixCSC{<:Real,Int64}#num_reactants*num_eqns
                          )::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1}
                          #Array{Tuple{Num_(reac+prod),List_stoich,List_ind}}
    reactants_list=Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1}()
    for eqn_ind in 1:num_eqns
        indlist=zeros(Int64,16)
        stoichlist=zeros(Int8,15)
        reactant_inds=findall(!iszero,reactants_mtx[:,eqn_ind])
        num_stoichs=length(reactant_inds)
        @assert(num_stoichs<=15)#or it would break the static Array
        for i in 1:num_stoichs
            indlist[i]=reactant_inds[i]
            stoichlist[i]=reactants_mtx[reactant_inds[i],eqn_ind]
        end
        indvec=SVector{16,Int64}(indlist)
        stoichvec=SVector{15,Int8}(stoichlist)
        push!(reactants_list,(num_stoichs,stoichvec,indvec))
    end
    @assert(length(reactants_list)==num_eqns)
    return reactants_list
end