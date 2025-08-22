using LinearAlgebra
function stab_func_bi(f_DDE, f_tau, pars, xe::Vector, nd, N; diff=1, hopf=0) #interval, N; diff=1)#(f_DDE, f_tau, pars, xe, nd, tj, ti,N; diff=1_) #tj are the interpolation point
    ##COULD ADD INTERVAL ARGUMENT BUT WANT TO ENSURE THEY DO IT OVER [-TAU_MAX , 0]
    #inputs:
    #xe is the equilibrium point you want to find the stability of 
    #pars are the parameters - ensure the parameters match the ones that given you equilibrium xe
    #interval is the interval over which you want to find the interpolation points over (e.g. [-tau_max, 0]) #may not be needed if you can work it out in function
    #N is the number to create N+1 nodes for interpolation

    m=length(xe)
    Id=Matrix{Float64}(I,m,m)
    xevec=[xe for _ in 1:nd+1]
    tau=f_tau(xevec, pars)
    taumax=findmax(tau)[1] #returns highest delay
    ti=vcat(0, -tau) #values to evaluate at [0,-tau1,-tau2,...,-tau_nd]
    nti=length(ti) #number of evaluation points
    #tj_int=[-taumax,0.0] #interval of tj values

    tj=reverse((-taumax/2)*(cos.(pi*(0:N)'/N)[:].+1)) #creates tj (interpoltaion) values in form of chebyshev points of 2nd kind over interval [-tau_max, 0]
    ljvals=j_eval(tj, ti) #finds l_j(0), l_j(-tau1),...,l_j(-tau_nd) for j=1 to N+1 (or j=0 to N in literature)
    D1=j_eval(tj,ti,diff=1) #finds first derivative matrix for interpolation points tj
    #stabmat= fill(NaN, m*(N+1), m*(N+1))#creates blank matrix for the spectral differentiation matrix A_N
    stabmat= fill(0.0, m*(N+1), m*(N+1))#creates blank matrix for the spectral differentiation matrix A_N

    mDmat=fill(NaN, m*N, m*(N+1)) #creates blank matrix for finding d_ij for i=2,...,N+1, j=1,...,N+1 (bottom rows of stabmat)

    for i in 2:N+1
        for j in 1:N+1
            mDmat[1+m*(i-2):m*(i-1),1+m*(j-1):m*j]=Id*D1[i,j]
        end 
    end 

    stabmat[m+1:end,1:end]=mDmat #fills in d_ij values of A_N matrix (stability matrix)

    function df(i,x) #function finds the partial derivative matrices (Jacobians)
        params=deepcopy(pars)
        J=f_deriv(f_DDE,x,params,nd,nx=i) #finds Jacobian for derivatives wrt x(t) for i=1, x(t-τ_1) for i=2, etc.
        return J
    end

    A=fill(NaN,m,m,nd+1) #creates array

    for i in 1:nd+1
        A[:,:,i]=df(i,xe) #finds A_0, A_1,...,A_nd
    end 
    
    for j in 1:N+1
        for k in 1:nti
            stabmat[1:m,1+m*(j-1):m*j]+=A[:,:,k]*(ljvals[k,j]*Id)
        end 
    end 
    #return A,taumax
    #return mDmat
    #return stabmat
    ev=eigen(stabmat)
    sm_eigvals=ev.values
    lambda_r_indx=findmax(real(sm_eigvals))[2] #finds index of rightmost eigenvalue
    lambda_r=sm_eigvals[lambda_r_indx] #returns rightmost eigenvalue

    stab=NaN

    if (real(lambda_r)>0.0) #equilibrium is unstable if real part of rightmost eigenvalue greater than 0
        stab=0
    else 
        stab=1 #equilibrium is stable if real part of rightmost eigenvalue less than 0
    end

    if hopf==0
        return stab, lambda_r
    else
        imin=argmin(abs.(real.(sm_eigvals))) #finds index of eigenvalue that is closest to being purely imaginary
    
        evecs=ev.vectors[:,imin]
        v0=evecs[1:m] #finds first eigenvector of eigenvalue closest to being purely imaginary

        #Seperates v0 into its real and imaginary parts
        vrini= real(v0) 
        viini=-imag(v0)

        #Below rescales vrini and viini so the norm is 1
        nv=norm(vcat(vrini,viini))
        vrini=vrini/nv
        viini=viini/nv

        omini=abs(imag(sm_eigvals[imin]))
        return stab, vrini, viini, omini #, lambda_r
    end

end 