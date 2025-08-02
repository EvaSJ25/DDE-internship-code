using LinearAlgebra #needed for I in identity matrix
function create_hopffunc(f_DDE,f_tau, pars, x0,p0,par_indx,nd;m=100) 
    include("f_deriv.jl") 
    include("stab_func.jl")
    n=length(x0) #number of states of x:x_1,...,x_n
    uvec1=[x0 for _ in 1:nd+1]
    params=deepcopy(pars)
    params[par_indx]=p0
    Id=Matrix{Float64}(I,n,n)

    function df(s,x,p)
        params=deepcopy(pars)
        params[par_indx]=p
        J=f_deriv(f_DDE,x,params,nd,nx=s)
        return J
    end

    sfunc=stab_func(f_DDE,f_tau,x0,p0,params,par_indx,nd,hopf=1) #want hopf!=0 here as we need the omega value
    vrini=sfunc[2]
    viini=sfunc[3]
    omini=sfunc[4]

    y0=vcat(x0,vrini,viini,omini,p0)

    function fhopf(y)
        params=deepcopy(pars)
        u,vr,vi,om,p=y[1:n], y[n+1:2*n], y[2*n+1:3*n],y[3*n+1],y[3*n+2] #u,vr,vi,om,p=y[1:2], y[3:4], y[5:6],y[7],y[8]
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
            A_mat+=A[:,:,j]*exp(-(om*im)*tau[j-1])
        end

        J=(om*im).*Id-A_mat #characteristic (equation) matrix J=(Δ(λ)) 
        rdf=J*v#charactersitic equation matrix times v (rdf=Δ(λ)*v) # v corresponds to v in system solution u(t)=ve^λt
        real_rdf=real(rdf)
        imag_rdf=imag(rdf)

        return vcat(rf,real_rdf,imag_rdf,vr'*vr+vi'*vi-1,viini'*vr-vrini'*vi)
    end 
    return y0, fhopf
end 