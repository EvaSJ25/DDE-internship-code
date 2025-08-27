using LinearAlgebra
function stab_func_DDE(A,taus, N; eigvecs=0)
    ##inputs: 
    #A is a vector of matrices, where A_0,...,A_d are the partial derivative matrices or the matrices in a linear DDE system
    #τ are the delays (in vector form)
    #N are the number of nodes for interpolation (so you get N+1 nodes)
    #eigvecs tells whether or not to output the eigenvectors (0 means no eigenvectors outputted, any other value means they are)

    ##outputs:
    #stab is the stability of the (equilibrium) point (already incorporated in the A matrices)
    #sm_eigvals are the eigenvalues of the point
    #if asked: sm_eigvecs are the eigenvectors (in matrix form) of the point - ith column of matrix contain the eigenvectors of the ith eigenvalue

    m=map(size,A)[1][1] #finds the number of states in the system 
    nA=length(A) #number of matrices 
    Id=Matrix{Float64}(I,m,m) #creates mxm identity matrix
    taumax=findmax(taus)[1] #returns highest delay
    ti=vcat(0, -taus) #values to evaluate at [0,-tau1,-tau2,...,-tau_nd]
    nti=length(ti) #number of evaluation points


    tj=reverse((-taumax/2)*(cos.(pi*(0:N)'/N)[:].+1)) #creates tj (interpoltaion) values in form of chebyshev points of 2nd kind over interval [-tau_max, 0]
    ljvals=j_eval(tj, ti) #finds l_j(0), l_j(-tau1),...,l_j(-tau_nd) for j=1 to N+1 (or j=0 to N in literature)
    D1=j_diff(tj) #finds first derivative matrix for interpolation points tj
    stabmat= fill(0.0, m*(N+1), m*(N+1))#creates blank matrix for the spectral differentiation matrix A_N
 
    mDmat=fill(NaN, m*N, m*(N+1)) #creates blank matrix for finding d_ij for i=2,...,N+1, j=1,...,N+1 (bottom rows of stabmat) - adjusted layout from (Breda et al 2009)

    for i in 2:N+1
        for j in 1:N+1
            mDmat[1+m*(i-2):m*(i-1),1+m*(j-1):m*j]=Id*D1[i,j] #finds d_ij values
        end 
    end 

    stabmat[m+1:end,1:end]=mDmat #fills in d_ij values of A_N matrix (stability matrix)

    for j in 1:N+1
        for k in 1:nA
            stabmat[1:m,1+m*(j-1):m*j]+=A[k]*(ljvals[k,j]*Id) #finds top row(s) of A_N - a_j for j=1,...,N+1 (adjusted from layout in (Breda et al 2009))
        end 
    end 

    ev=eigen(stabmat)
    sm_eigvals=ev.values #stability matrix's eigenvalues
 
    lambda_r_indx=findmax(real(sm_eigvals))[2] #finds index of rightmost eigenvalue
    lambda_r=sm_eigvals[lambda_r_indx] #returns rightmost eigenvalue

    stab=NaN #creates blank stability indicator

    if (real(lambda_r)>0.0) #equilibrium is unstable if real part of rightmost eigenvalue greater than 0
        stab=0
    else 
        stab=1 #equilibrium is stable if real part of rightmost eigenvalue less than 0
    end
    
    if eigvecs!=0 #if function asked for eigenvectors
        sm_eigvecs=ev.vectors #stability matrix's eigenvectors
        return stab, sm_eigvals, sm_eigvecs #returns stability, eigenvalues and eigenvectors
    else 
        return stab, sm_eigvals #only returns stability and eigenvalues
    end

    #Reference:
    #D. Breda, S. Maset, R. Vermiglio. *TACE-DDE: a Tool for Robust Analysis and Characteristic Equations for Delay Differential Equations*, volume 388, pages 145-155, 2009
    
end 