function stab_func(f_DDE,f_tau,x0::Vector,p0::Vector,pars,par_indx::Vector,nd;doprint=1,hopf=0,h=1e-6,m=100)#(f_DDE,f_tau,x0,p0,pars,par_indx,nd;doprint=1,hopf=0,h=1e-6,m=100)
    ##inputs:
    #f_DDE is the DDE system
    #f_tau is the function for the delay
    #x0 is the equilibirium point you're finding the stability of
    #p0 is the (varied) parameter value at which this equilibrium point occurs at
    #pars= include constant tau value
    #m in the number of steps you want (to discretise over)
    #pars_indx is the parameter that was varied to get xlist, plist (using track_curve function)
    #nd is the number of delays
    #doprint=1 means that the eigvalues and the lowest eigenvector is given, doprint=0 doesn't return these
    #hopf=1 outputs the estimated ω value if the user wants a hopf bifurcation. hopf=0 doesn't return ω value

    ##outputs:
    #stab=the stability of the equilibrium point given
    #eigvalsJ=the eigenvalues of the equilibrium point
    #eigvecs=the eigenvectors of the equilibrium point
    ##if hopf!=0 then the function can also output
    #vrini=initial guess of real part of first eigenvector of eigenvalue closest to being purely imaginary 
    #viini=initial guess of imaginary part of first eigenvector of eigenvalue closest to being purely imaginary 
    #omini=initial guess of ω

    n=length(x0) #number of states
    uvec1=[x0 for _ in 1:nd+1] #repeats the equilibrium point for all delayed states
    params=deepcopy(pars)
    params[par_indx]=p0 
    Id=Matrix{Float64}(I,n,n) #creates identity matrix

    l=n*(1+nd*m) #number of rows (=columns) of large matrix
    
    tau=f_tau(uvec1,params) #finds tau(s) of system (so they can be used in the large matrix calculation)
    stab_mat=fill(0.0,l,l) #we assume that m_1=m_2=m (e.g. for nd=2)

    #-m/tau1 on the diagonal 
    stab_mat[diagind(stab_mat)[1:n*(1+m)]].=-(m/tau[1])

    #puts -m/tau2,...,-m/tau_nd on leading diagonal in the right places for nd delays
    for i in 2:nd 
        stab_mat[diagind(stab_mat)[n*(1+(i-1)*m)+1:n*(1+i*m)]].=-(m/tau[i])
    end 

    #puts m/tau_1,...,m/tau_nd on the right spots on the off diagonal
    for k in 1:nd
        for j in n*((k-1)*m+1)+1:n*(k*m+1)
            stab_mat[j,j-n]=m/tau[k]
        end 
    end
    
    #For m/tau_k, on the n((k-1)m+1)+1:n((k-1)m+1)+n rows the (m/tau_k)*Id is in the first columns and the off diagonal entries are 0
    for i in 2:nd
        stab_mat[n*((i-1)*m+1)+1:n*((i-1)*m+1)+n,1:n]=m/tau[i]*Id #in first columns m/tau_k is there (so v0 is in multiplication -see overleaf)
        stab_mat[n*((i-1)*m+1)+1:n*((i-1)*m+1)+n,n*((i-1)*m+1)+1-n:n*((i-1)*m+1)+n-n]=0*Id
    end 

    function df(i,x,p) #function finds the partial derivative matrices (Jacobians)
        params=deepcopy(pars)
        params[par_indx]=p
        J=f_deriv(f_DDE,x,params,nd,nx=i) #finds Jacobian for derivatives wrt x(t) for i=1, x(t-τ_1) for i=2, etc.
        return J
    end
    
    A=fill(NaN,n,n,nd+1)

    for i in 1:nd+1
        A[:,:,i]=df(i,x0,p0) #finds A_0, A_1,...,A_nd
    end 
    
    #adds A_0,A_1,...,A_nd matrices
    for j in 1:nd+1 #j in 1:nd-1
        stab_mat[1:n,(j-1)*m*n+1:n*((j-1)*m+1)]=A[:,:,j]
    end

    ev=eigen(stab_mat)
    eigvalsJ=ev.values #eigenvalues of stability matrix
    eigvecs=ev.vectors #eigenvectors of stability matrix
    
    if hopf!=0 #if you're looking at a Hopf bifurcation
        imin=argmin(abs.(real.(eigvalsJ))) #finds index of eigenvalue that is closest to being purely imaginary
    
        evecs=ev.vectors[:,imin]
        v0=evecs[1:n] #finds first eigenvector of eigenvalue closest to being purely imaginary

        #Seperates v0 into its real and imaginary parts
        vrini= real(v0) 
        viini=-imag(v0)

        #Below rescales vrini and viini so the norm is 1
        nv=norm(vcat(vrini,viini))
        vrini=vrini/nv
        viini=viini/nv

        omini=abs(imag(eigvalsJ[imin])) #initial guess for ω
    end

    lamb_dom1,=findmax(real(eigvalsJ))#finds the eigenvalue with the largest real part
    #if this lamb_dom1 is less than 0 then the equilibria is stable else it is unstable 
    
    if (lamb_dom1>0.0) #equilibrium is unstable
        stab=0
    else 
        stab=1 #equilibrium is stable
    end

    if doprint==1 #wants stability AND eigenvalues, eigenvectors (and more for Hopf)
        if hopf==0 #if not interested in hopf just return the stability and eigenvalues
            return stab, eigvalsJ,eigvecs
        else #return relevant hopf information
            return stab, vrini, viini, omini
        end
    else #if doprint is 0 then just return the stability of the point (no eigenvalues, no eigenvectors are returned)
        return stab
    end 
end 