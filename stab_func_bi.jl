function stab_func_bi(f_DDE, f_tau, pars, xe::Vector,p0::Vector, par_indx, nd, N; diff=1) #interval, N; diff=1)#(f_DDE, f_tau, pars, xe, nd, tj, ti,N; diff=1_) #tj are the interpolation point
    #inputs:
    #interval is the interval over which you want to find the interpolation points over (e.g. [-tau_max, 0]) #may not be needed if you can work it out in function
    #N is the number to create N+1 nodes for interpolation

    ndim=length(xe)
    xevec=[xe for _ in 1:nd+1]
    tau=f_tau(xevec, pars)
    taumax=findmax(tau)[1] #returns highest delay
    ti=vcat(0, -tau) #values to evaluate at [0,-tau1,-tau2,...,-tau_nd]
    tj_int=[-taumax,0.0] #interval of tj values

    function df(i,x,p) #function finds the partial derivative matrices (Jacobians)
        params=deepcopy(pars)
        params[par_indx]=p
        #@infiltrate
        J=f_deriv(f_DDE,x,params,nd,nx=i) #finds Jacobian for derivatives wrt x(t) for i=1, x(t-τ_1) for i=2, etc.
        return J
    end

    A=fill(NaN,ndim,ndim,nd+1) #creates array

    for i in 1:nd+1
        #@infiltrate
        A[:,:,i]=df(i,xe,p0) #finds A_0, A_1,...,A_nd
    end 
    
    return A,tau

end 