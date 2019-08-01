
function mk_reactants_list(num_reactants::Int,num_eqns::Int,
                           reactants_mtx::SparseMatrixCSC{Float64,Int64}#num_reactants*num_eqns
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

function generate_loss_gain(num_reactants::Int,num_eqns::Int,
                         stoich_mtx::SparseMatrixCSC{Float64,Int64})#num_reactants*num_eqns
    rate_prod_exprs=quote end#Array{Expr,1}()#num_eqns
    for eqn_ind in 1:num_eqns
        rate_expr=:(*(rate_values[$eqn_ind]))#assume there is at least one reactant
        reactant_inds=findall(!iszero,stoich_mtx[:,eqn_ind])
        @assert(length(reactant_inds)>0)
        for reactant_ind in reactant_inds
            stoich=Int(stoich_mtx[reactant_ind,eqn_ind])
            if stoich==1
                push!(rate_expr.args,:(reactants[$reactant_ind]))#r=[A]*[B]*...
            elseif stoich>=2
                push!(rate_expr.args,:(reactants[$reactant_ind]^$stoich))#r=[A]^a*[B]^b*...
            end
        end
        rate_prod=:(rate_prods[$eqn_ind]=$rate_expr)
        push!(rate_prod_exprs.args,rate_prod)
    end
    stoich_mtx_T=transpose(stoich_mtx)#num_eqns*num_reactants, column indexing is faster for CSC format
    loss_gain_exprs=quote end
    for reactant_ind in 1:num_reactants
        loss_gain_expr=:(dydt[$reactant_ind]=(+(0)))#(:(=), (:ref, :dydt, $reactant_ind), (:call, :+))
        eqn_inds=findall(!iszero,stoich_mtx_T[:,reactant_ind])
        @assert(length(eqn_inds)>0)
        for eqn_ind in eqn_inds
            coef=Int((-1)*stoich_mtx_T[eqn_ind,reactant_ind])#negative for reactants, positive for products
            #rate_expr=rate_exprs[eqn_ind]
            if coef==1
                push!(loss_gain_expr.args[2].args,:(rate_prods[$eqn_ind]))
            else
                push!(loss_gain_expr.args[2].args,:($coef*rate_prods[$eqn_ind]))
            end
        end
        push!(loss_gain_exprs.args,loss_gain_expr)
    end
    return quote
        function loss_gain_static!(num_reactants::Int,num_eqns::Int,
                                   reactants::Array{Float64,1},#num_reactants
                                   rate_values::Array{Float64,1},#num_eqns
                                   rate_prods::Array{Float64,1},#num_eqns k*[A]^a*[B]^b
                                   dydt::Array{Float64,1}#num_reactants
                                  )
                $rate_prod_exprs
                $loss_gain_exprs
            return dydt
        end
    end
end

function constant_folding!(fun_expr::Expr,constant_list::Dict,rate_values::Array{<:Real,1})
    for (key, value) in constant_list
        eval(:($key=$value))
    end
        
    function try_eval(expr::Expr)::Expr
        #expr :(x=1+1) => :(x=2), :(x=a+b) => :(x=a+b), :(x=1+y),y=2 =>:(x=3)
        try
            ans=eval(expr.args[2])
            var=expr.args[1]
            return :($var=$ans)
        catch
            return expr
        end
    end
        
    ops=fun_expr.args[2].args[2].args
    for i in 1:length(ops)
        op=ops[i]
        if typeof(op)==LineNumberNode
            continue
        end
        if op.head==:(=)
            expr=try_eval(op)
            if typeof(expr.args[2])<:Number
                eval(expr) #bring var into scope
                ops[i]=expr
            end
        end
        if op.head==:block
            rate_val_ops=op.args
            for j in 1:length(rate_val_ops)
                r_op=rate_val_ops[j]
                if typeof(r_op)==LineNumberNode
                    continue
                end
                if r_op.head==:(=)
                    expr=try_eval(r_op)
                    if typeof(expr.args[2])<:Number
                        rate_val_ops[j]=expr #extract constant values
                        #eval(expr) #like rate_values[5]=20. update rate_values if constant
                        rate_values[expr.args[1].args[2]]=expr.args[2]#expr :(ratevalues[5]=4)<=>(:(=), (:ref, :ratevalues, 5), 4)
                    end
                end
            end
        end        
    end
    return fun_expr
end

function extract_constants!(fun_expr::Expr)::Expr
    ops=fun_expr.args[2].args[2].args
    for op in ops
        if typeof(op)==LineNumberNode
            continue
        end
        if op.head==:block
            filter!(l->!(typeof(l)==LineNumberNode),op.args)
            filter!(l->!(typeof(l.args[2])<:Number),op.args)
        end
    end
    return fun_expr
end