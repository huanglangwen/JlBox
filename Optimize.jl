
function constant_folding!(fun_expr::Expr,constant_list::Dict,rate_values::Array{Float64,1})
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
        if op.head==:block
            filter!(l->!(typeof(l.args[2])<:Number),op.args)
        end
    end
    return fun_expr
end