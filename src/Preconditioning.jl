using IncompleteLU

struct TridiagonalLUCache
    LU::Tridiagonal# L: LU.dl & 1, U: lu.d & lu.du
    cache::Array{Float64,1} 
end

function factorize!(lu::TridiagonalLUCache,trimat::Tridiagonal)
    lu.LU.d[1] = trimat.d[1]
    n, _ = size(trimat)
    for i = 2:n
        dll = trimat.dl[i-1]/lu.LU.d[i-1]
        lu.LU.d[i] = trimat.d[i] - dll*trimat.du[i-1]
        lu.LU.dl[i-1] = dll
    end
    lu.LU.du .= trimat.du
end

function LinearAlgebra.ldiv!(x,A::TridiagonalLUCache,y)
    A.cache[1] = y[1]
    n = length(y)
    for i = 2:n
        A.cache[i] = y[i] - A.LU.dl[i-1]*A.cache[i-1]
    end
    x[n] = A.cache[n]/A.LU.d[n]
    for i = n-1:-1:1
        x[i] = (A.cache[i] - A.LU.du[i]*x[i+1])/A.LU.d[i]
    end
end

function sp2tridiag(spmat, gamma)
    n, _ = size(spmat)
    dl = zeros(n)
    du = zeros(n)
    d = zeros(n)
    for j = 1:n
        dl[j] = -gamma*spmat[min(j+1,n),j]
        d[j] = 1-gamma*spmat[j,j]
        du[j] = -gamma*spmat[max(j-1,1),j]
    end
    return Tridiagonal(dl[1:end-1],d,du[2:end])
end

function sp2tridiag!(trimat, spmat, gamma)
    n, _ = size(spmat)
    for j = 1:n
        if (j+1 <= n) trimat.dl[j] = -gamma*spmat[j+1,j] end
        trimat.d[j] = 1-gamma*spmat[j,j]
        if (j >= 2) trimat.du[j-1] = -gamma*spmat[j-1,j] end
    end
end

function J2M!(spmat::SparseMatrixCSC, gamma)#M=I-gamma*J
    n, _ = size(spmat)
    rows = rowvals(spmat)
    vals = nonzeros(spmat)
    for j = 1:n
        for i in nzrange(spmat,j)
            id = rows[i] == j ? 1 : 0
            vals[i] = id - gamma*vals[i]
        end
    end
end

function J2M!(trimat::Tridiagonal, gamma)#M=I-gamma*J
    trimat.dl .*= -gamma
    trimat.du .*= -gamma
    trimat.d .= 1.0 .- gamma.*trimat.d
end

function default_prec()
    return (z,r,p,t,y,fy,gamma,delta,lr)->LinearAlgebra.ldiv!(z,p["preconditioner"],r)
end

function default_psetup(jacmethod1::String, jacmethod2::String, initsteps::Int = 20)
    jac1 = select_jacobian(jacmethod1, nothing)
    jac2 = select_jacobian(jacmethod2, nothing)
    return function (p,t,u,du,jok,jcurPtr,gamma)#factorize(I-gamma*J)
        if p["Current_iter"] < initsteps
            spmat = p["sparsity"]
            jac1(spmat,u,p,t)
            J2M!(spmat,gamma)
            p["preconditioner"]=IncompleteLU.ilu(spmat;τ=1e-5)
            @debug "PSETUP CALLED AS INIT"
        elseif !(p["preconditioner"] isa TridiagonalLUCache)
            n,_ = size(p["sparsity"])
            jac2(p["sparsity"],u,p,t)
            p["trimat"]=sp2tridiag(p["sparsity"],gamma)
            p["preconditioner"]=TridiagonalLUCache(p["trimat"],zeros(n))
            factorize!(p["preconditioner"],p["trimat"])
        else
            jac2(p["sparsity"],u,p,t)
            sp2tridiag!(p["trimat"],p["sparsity"],gamma)
            factorize!(p["preconditioner"],p["trimat"])
            @debug "PSETUP CALLED"
        end
        jcurPtr[]=false
    end
end

function default_psetup(jacmethod, initsteps::Int = 20)
    return default_psetup(jacmethod, jacmethod, initsteps)
end

function default_psetup_gas()
    return function (p,t,u,du,jok,jcurPtr,gamma)#factorize(I-gamma*J)
        if !(p["sparsity"] isa Tridiagonal)
            n,_ = size(p["sparsity"])
            p["sparsity"] = sp2tridiag(p["sparsity"],gamma)
            p["preconditioner"] = TridiagonalLUCache(p["sparsity"],zeros(n))
        end
        trimat = p["sparsity"]
        gas_jac_banded!(trimat,u,p,t)
        J2M!(trimat,gamma)
        factorize!(p["preconditioner"],trimat)
        @debug "PSETUP CALLED"
        jcurPtr[]=false
    end
end

function gas_jac_banded!(jac_mtx,reactants::Array{<:Real,1},p::Dict,t::Real)
    rate_values,stoich_mtx,stoich_list,reactants_list,num_eqns,num_reactants=
        [p[ind] for ind in 
            ["rate_values","stoich_mtx","stoich_list","reactants_list",
             "num_eqns","num_reactants"]
        ]
    #Probably have to re-eval rate_values again
    loss_gain_jac_banded!(-1,1,num_reactants,num_eqns,reactants,stoich_mtx,stoich_list,reactants_list,rate_values,jac_mtx)
end

function loss_gain_jac_banded!(
                       lower::Int,upper::Int,#Lower & upper band of jacobian matrix
                       num_reactants::Int,num_eqns::Int,
                       reactants::Array{<:Real,1},#num_reactants
                       stoich_mtx::SparseMatrixCSC{Float64,Int64},#num_reactants*num_eqns
                       stoich_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, both reac and prod
                       reactants_list::Array{Tuple{Int8,SVector{15,Int8},SVector{16,Int64}},1},#num_eqns, only reac
                       rate_values::Array{<:Real,1},#num_eqns
                       lossgain_jac_mtx#::Array{Float64,2}#SparseMatrixCSC{Float64,Int64},#num_output(dydt)*num_input(y)
                       )
    #lossgain_jac_mtx=spzeros(num_reactants,num_reactants)#num_output(dydt)*num_input(y)
    fill!(lossgain_jac_mtx, 0.0)
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
                band=reactant_y_ind - reactant_ind
                if band >= lower && band <= upper
                    lossgain_jac_mtx[reactant_ind,reactant_y_ind]+=stoich_mtx[reactant_ind,eqn_ind]*prod*(-1)
                end
            end
        end
    end
    nothing
end