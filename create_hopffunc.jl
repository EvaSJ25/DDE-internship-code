using LinearAlgebra #needed for I in identity matrix
function create_hopffunc(f_DDE,f_tau, pars, x0,p0::Vector,par_indx::Vector,nd;m=100) 
    ##inputs:
    #f_DDE is the system function(s)
    #f_tau is the delay equation/function
    #pars are the parameters values of the system
    #x0 is the initial guess of the (equilibrium) state the Hopf bifurcation
    #p0 is the initial parameter guess for the Hopf bifurcation (note this always a vector (even if only varying one parameter))
    #par_indx is the parameter(s) you're varying (again it should be given as a vector even if only of length 1)
    #nd is the number of delays
    #m is the number of discretised steps to be used stab_func
    
    ##outputs:
    #y0 are the initial guesses of x, vr, vi, om, p
    #Note:where x is are the states, vr and vi are the real and imaginary parts of the (first) eigenvector for eigenvalue om (ω) and p is the initial guess of Hopf parameter p
    #fhopf is the function that finds the Hopf and is in a form that can be used in Newton (for finding exact parameter values) or Track curve (to continue to 2-parameter plane)

    include("f_deriv.jl") #used to find the state partial derivatives matrices
    include("stab_func.jl") #used to find the eigenvalues and eigenvectors of initial Hopf guess
    n=length(x0) #number of states of x:x_1,...,x_n
    uvec1=[x0 for _ in 1:nd+1] #makes vector of vectors (repeats equilibrium for x(t)=x(t-τ1)=..=x(t-τnd)
    params=deepcopy(pars)
    params[par_indx]=p0 #sets varied parameter to parameter guess p0

    Id=Matrix{Float64}(I,n,n) #sets up identity matrix

    #The below function is for finding the partial derivative matrices (A_0,A_1,...,A_nd)
    function df(i,x,p) #note i is the desired state derivative (i=1 is wrt x(t), i=2 is wrt x(t-tau_1), etc.)
        params=deepcopy(pars)
        params[par_indx]=p
        J=f_deriv(f_DDE,x,params,nd,nx=i)
        return J
    end

    sfunc=stab_func(f_DDE,f_tau,x0,p0,params,par_indx,nd,hopf=1) #want hopf!=0 here as we need the omega value
    vrini=sfunc[2] #real part of (first) eigenvector for eigenvalue closest to being purely imaginary
    viini=sfunc[3] #imaginary part of (first) eigenvector for eigenvalue closest to being purely imaginary
    omini=sfunc[4] #initial ω guess (the eigenvalue closest to being purely imaginary)
    

    y0=vcat(x0,vrini,viini,omini,p0) #combines initial guess of relevant Hopf information

    #Below defines the function for finding the Hopf (involving the characteristic equation, etc.)
    function fhopf(y) 
        params=deepcopy(pars)
        u,vr,vi,om,p=y[1:n], y[n+1:2*n], y[2*n+1:3*n],y[3*n+1],y[3*n+2:end] 
        params[par_indx]=p

        uvec=[u for _ in 1:nd+1]
        v=vr+vi*im
        A=fill(NaN,n,n,nd+1)

        for i in 1:nd+1
            A[:,:,i]=df(i,u,p)
        end 

        rf=f_DDE(uvec,params) # still need f(x*,...x*,params)=0
        A_mat=A[:,:,1] #starts with A_0

        tau=f_tau(uvec,params)
        
        #Below creates the sum A_0 + A_1*exp(-lam*tau_1)+...+ A_nd*(-lam*tau_nd)
        for j in 2:nd+1
            A_mat+=A[:,:,j]*exp(-(om*im)*tau[j-1])
        end

        J=(om*im).*Id-A_mat #characteristic (equation) matrix J=(Δ(λ)) 
        rdf=J*v #charactersitic equation matrix times v (rdf=Δ(λ)*v) # v corresponds to v in system solution u(t)=ve^λt
        real_rdf=real(rdf)
        imag_rdf=imag(rdf)

        return vcat(rf,real_rdf,imag_rdf,vr'*vr+vi'*vi-1,viini'*vr-vrini'*vi)
    end 
    return y0, fhopf
end 