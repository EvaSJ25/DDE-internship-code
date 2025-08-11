using LinearAlgebra
function create_foldfunc(f_DDE, f_tau,pars,x0,p0,par_indx,nd;m=100)
    #find eigenvalue that is closest to being purely 0
    include("stab_func.jl")
    include("f_deriv.jl")

    n=length(x0) #number of states of x:x_1,...,x_n
    np=length(par_indx)
    uvec1=[x0 for _ in 1:nd+1]
    params=deepcopy(pars)
    params[par_indx]=p0
    #Id=Matrix{Float64}(I,n,n) #not needed as lambda is 0 for a fold 

    function df(s,x,p) #finds the partial derivative matrices
        params=deepcopy(pars)
        params[par_indx]=p
        J=f_deriv(f_DDE,x,params,nd,nx=s)
        return J
    end

    stabfunc=stab_func(f_DDE, f_tau,x0,p0,params,par_indx,nd) #returns stability, eigenvalues and eigenvectors for point x0
    eigvals=stabfunc[2] #eigenvalues 
    eigvecs=stabfunc[3] #eigenvectors

    fold_indx=argmin(abs.(eigvals))#finds index of the eigenvalue that is closest to being purely 0
    #lamini=real(eigvals[fold_indx]) #eigenvalue closest to being 0
    evecs=eigvecs[:,fold_indx]
    v0=evecs[1:n]

    vrini=real(v0)
    viini=-imag(v0)
    nv=norm(vcat(vrini,viini))
    vrini=vrini/nv
    viini=viini/nv

    #y0=vcat(x0,vrini,viini,lamini,p0)
    #@infiltrate
    y0=vcat(x0,vrini,viini,p0)
    function ffold(y)
        params=deepcopy(pars)
        if length(par_indx)==1
            #u,vr,vi,lam,p=y[1:n], y[n+1:2*n], y[2*n+1:3*n],y[3*n+1],y[3*n+2] #u,vr,vi,om,p=y[1:2], y[3:4], y[5:6],y[7],y[8]
            u,vr,vi,p=y[1:n], y[n+1:2*n], y[2*n+1:3*n],y[3*n+1] #u,vr,vi,om,p=y[1:2], y[3:4], y[5:6],y[7],y[8]

        else
            #u,vr,vi,lam,p=y[1:n], y[n+1:2*n], y[2*n+1:3*n],y[3*n+1],y[3*n+2:end] #u,vr,vi,om,p=y[1:2], y[3:4], y[5:6],y[7],y[8]
            u,vr,vi,lam,p=y[1:n], y[n+1:2*n], y[2*n+1:3*n],y[3*n+1:end] #u,vr,vi,om,p=y[1:2], y[3:4], y[5:6],y[7],y[8]
        end
        #u,vr,vi,p=y[1:n],y[n+1:2*n],y[2*n+1:3*n],y[3*n+1]
        #u,vr,vi,p=y[1:n], y[n+1:2*n], y[2*n+1:3*n],y[3*n+1],y[3*n+2] #u,vr,vi,om,p=y[1:2], y[3:4], y[5:6],y[7],y[8]
        params[par_indx]=p
        uvec=[u for _ in 1:nd+1]
        v=vr+vi*im

        A=fill(NaN,n,n,nd+1)
        for i in 1:nd+1
            A[:,:,i]=df(i,u,p)
        end 

        rf=f_DDE(uvec,params)
        A_mat=A[:,:,1] #starts with A_0
        #Below creates the sum A_0 + A_1*exp(-lam*tau_1)+...+ A_nd*(-lam*tau_nd)
        for j in 2:nd+1
            A_mat+=A[:,:,j]
        end

        J=-A_mat #characteristic (equation) matrix J=(Δ(λ)) for fold bifurcation (where λ=0)
        #J=lam.*Id-A_mat
        rdf=J*v#charactersitic equation matrix times v (rdf=Δ(λ)*v) # v corresponds to v in system solution u(t)=ve^λt
        real_rdf=real(rdf)
        imag_rdf=imag(rdf)

        return vcat(rf,real_rdf,imag_rdf,vr'*vr+vi'*vi-1)
    end 

    return y0, ffold

end 