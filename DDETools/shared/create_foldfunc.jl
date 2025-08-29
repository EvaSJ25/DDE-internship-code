using LinearAlgebra
function create_foldfunc(f_DDE, f_tau,pars,x0,p0::Vector,par_indx::Vector,nd;m=100) #helps to initialise guess and function to find eigenvalue that is closest to being purely 0
    ##inputs:
    #f_DDE is the system function(s)
    #f_tau is the delay equation/function
    #pars are the parameters values of the system
    #x0 is the initial guess of the (equilibrium) state the fold bifurcation
    #p0 is the initial parameter guess for the fold bifurcation (note this should always be given as a vector)
    #par_indx is the index of the parameter(s) you're varying (again it should be given as a vector even if only of length 1)
    #nd is the number of delays
    #m is the number of discretised steps to be used in stab_func_matrix
    
    ##outputs:
    #y0 are the initial guesses of x, vrini, viini, p
    #Note:where x is are the states, vrini and viini are the real and imaginary parts of the (first) eigenvector for eigenvalue closest to being 0 and p is the initial guess of fold parameter p
    #ffold is the function that finds the fold bifurcation


    n=length(x0) #number of states of x (i.e. x_1,...,x_n)
    uvec1=[x0 for _ in 1:nd+1] #makes vector of vectors (repeats equilibrium for x(t)=x(t-τ1)=..=x(t-τnd))
    params=deepcopy(pars)
    params[par_indx]=p0 #sets varied parameter to parameter guess p0

    #The below function is for finding the partial derivative matrices (A_0,A_1,...,A_nd)
    function df(i,x,p) #note i is the desired state derivative (i=1 is wrt x(t), i=2 is wrt x(t-tau_1), etc.)
        params=deepcopy(pars)
        params[par_indx]=p
        J=f_deriv(f_DDE,x,params,nd,nx=i)
        return J
    end

    stabfunc=stab_func_matrix(f_DDE, f_tau,x0,p0,params,par_indx,nd) #returns stability, eigenvalues and eigenvectors for point x0
    eigvals=stabfunc[2] #eigenvalues 
    eigvecs=stabfunc[3] #eigenvectors

    fold_indx=argmin(abs.(eigvals))#finds index of the eigenvalue that is closest to being purely 0
    evecs=eigvecs[:,fold_indx] #eigenvectors of eigenvalue closest to 0
    v0=evecs[1:n] #first eigenvector of eigenvalue closest to 0

    #Seperates v0 into its real and imaginary parts
    vrini=real(v0)
    viini=-imag(v0)

    #Below rescales vrini and viini so that their norm equals 1
    nv=norm(vcat(vrini,viini))
    vrini=vrini/nv
    viini=viini/nv

    y0=vcat(x0,vrini,viini,p0) #combines initial guess of relevant fold information

    #Below defines the function for finding the fold (involving the characteristic equation, etc.)
    function ffold(y)
        params=deepcopy(pars)
        u,vr,vi,p=y[1:n], y[n+1:2*n], y[2*n+1:3*n],y[3*n+1:end] 
        params[par_indx]=p
        uvec=[u for _ in 1:nd+1]
        v=vr+vi*im #combines the real and imaginary parts of the eigenvector

        A=fill(NaN,n,n,nd+1)
        for i in 1:nd+1
            A[:,:,i]=df(i,u,p) #Finds A_0,A_1,...,A_nd matrices
        end 

        rf=f_DDE(uvec,params) #still need f(x*,...x*,params)=0 condition for a fold bifurcation

        A_mat=A[:,:,1] #starts with A_0
        #Below creates the sum A_0 + A_1*exp(-λ*τ_1)+...+ A_nd*(-λ*τ_nd)=A_0 + A_1+..._A_nd (as λ=0)
        for j in 2:nd+1
            A_mat+=A[:,:,j]
        end

        J=-A_mat #characteristic (equation) matrix J=(Δ(λ)) for fold bifurcation (where λ=0)
        rdf=J*v#charactersitic equation matrix times v (rdf=Δ(λ)*v) # v corresponds to v in system solution u(t)=ve^λt
        real_rdf=real(rdf)
        imag_rdf=imag(rdf)

        return vcat(rf,real_rdf,imag_rdf,vr'*vr+vi'*vi-1) #these conditions all need to equal 0 to find the true fold bifurcation
    end 

    return y0, ffold

end 