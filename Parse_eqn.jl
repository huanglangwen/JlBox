using SparseArrays

function stoich_str2int(stoich::String)::Int
    if stoich==""
        return 1
    else
        return parse(Int,stoich)
    end
end

function parse_reactants(file::String)
    reactants_dict=Dict{String,Set{Tuple{Int,Int}}}() # Dict{Reactants,Set{Tuple{LineNum,Stoich}}}
    reactants2ind=Dict{String,Int}()
    ind2reactants=Dict{Int,String}()
    reactant_step=0
    rate_values=zeros(0)
    #file="MCM_test.eqn.txt"
    num_eqns=0
    open(file,"r") do f
        #global num_eqns
        ind=0
        for line in eachline(f)
            ind +=1
            num_eqns=max(ind,num_eqns)
            line_parts=split(strip(line),",")
            push!(rate_values,parse(Float64,line_parts[2]))
            eqn=line_parts[1]
            lhs_rhs=[strip(i) for i in split(eqn,"->")]
            lhs=[strip(i) for i in split(lhs_rhs[1],"+")]
            rhs=[strip(i) for i in split(lhs_rhs[2],"+")]
            reac_regexp=r"([0-9]*)([A-Z][0-9A-Z]*)"
            if lhs[1]==""#in some eqns, O2 and H2O is ignored, like O+O3=
                reacs=[]
            else
                reacs=[Tuple{String,String}(match(reac_regexp,i).captures) for i in lhs]
            end
            if rhs[1]==""
                prods=[]
            else
                prods=[Tuple{String,String}(match(reac_regexp,i).captures) for i in rhs]
            end
            for (stoich,reac) in reacs
                num_stoich=stoich_str2int(stoich)
                if !haskey(reactants_dict,reac)
                    reactants_dict[reac]=Set()
                    reactant_step+=1
                    reactants2ind[reac]=reactant_step
                    ind2reactants[reactant_step]=reac
                end
                push!(reactants_dict[reac],(ind,num_stoich))
            end
            for (stoich,reac) in prods
                num_stoich=stoich_str2int(stoich)*(-1)# because it is product
                if !haskey(reactants_dict,reac)
                    reactants_dict[reac]=Set()
                    reactant_step+=1
                    reactants2ind[reac]=reactant_step
                    ind2reactants[reactant_step]=reac
                end
                push!(reactants_dict[reac],(ind,num_stoich))
            end
            #println(reacs)
            #println(prods)
        end
    end

    num_reactants=length(reactants_dict)
    @assert(num_reactants==reactant_step)
    @assert(num_eqns==length(rate_values))
    #reactants_list=[ for i in 1:num_reactants]
    reactants_inds=1:num_reactants
    #reactants2ind=Dict(reactants_list[i]=>i for i in reactants_inds)

    stoich_mtx=spzeros(num_reactants,num_eqns)
    reactants_mtx=spzeros(num_reactants,num_eqns)
    for i in reactants_inds
        for (eqn_ind,stoich) in reactants_dict[ind2reactants[i]]
            stoich_mtx[i,eqn_ind]+=stoich#FOR DUPLICATED REACTANTS A+A->B+B
            if stoich>0
                reactants_mtx[i,eqn_ind]+=stoich #FOR CATALYSE A+B->A+C
            end
        end
    end
    dropzeros!(stoich_mtx)
    return (rate_values,stoich_mtx,reactants_mtx,num_eqns,num_reactants,reactants2ind)
end

function parse_init(file::String,reactants2ind::Dict,num_reactants::Integer)
    reactants_init=zeros(num_reactants)
    open(file,"r") do f
        for line in eachline(f)
            reac,val=split(line," = ")
            reactants_init[reactants2ind[reac]]=parse(Float64,val)
        end
    end
    reactants_init
end