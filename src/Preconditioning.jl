using IncompleteLU
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

function J2M(spmat, gamma)#M=I-gamma*J
    n, _ = size(spmat)
    rows = rowvals(spmat)
    vals = nonzeros(spmat)
    for j = 1:n
        for i in nzrange(spmat,j)
            id = rows[i] == j ? 1 : 0
            vals[i] = id - gamma*vals[i]
        end
    end
    return spmat
end

function default_prec()
    return (z,r,p,t,y,fy,gamma,delta,lr)->LinearAlgebra.ldiv!(z,p["preconditioner"],r)
end

function default_psetup(jacmethod1::String, jacmethod2::String, initsteps::Int = 20)
    jac1 = select_jacobian(jacmethod1, nothing)
    jac2 = select_jacobian(jacmethod2, nothing)
    return function (p,t,u,du,jok,jcurPtr,gamma)#factorize(I-gamma*J)
        if p["Current_iter"]< initsteps
            jac1(p["sparsity"],u,p,t)
            p["preconditioner"]=ilu(JlBox.J2M(p["sparsity"],gamma);τ=1e-5)
            @debug "PSETUP CALLED AS INIT"
        else
            jac2(p["sparsity"],u,p,t)
            p["preconditioner"]=factorize(JlBox.sp2tridiag(p["sparsity"],gamma))
            @debug "PSETUP CALLED"
        end
        jcurPtr[]=false
    end
end

function default_psetup(jacmethod, initsteps::Int = 20)
    return default_psetup(jacmethod, jacmethod, initsteps)
end