function f_deriv(f,u,pars,nd;nx=[],np=[],v=[],h=1e-5,k=1e-5) #u many columns
    ##inputs:
    #f is the system (equations)
    #u is the point the derivatives are being taken at
    #pars are the system parameters
    #nd is number of delays in the system
    #nx is the state derivative you want (empty=no state derivative wanted) -> nx=1 is derivtive wrt to state at time t, nx=2 is derivative wrt to state at t-τ_1, etc.
    #np is the index of the parameter derivative user wants (empty = no parameter derivative wanted)
    #h is the small step for finite difference 
    #k is the small steps for parameter (used only if nx and np both non-empty)
    
    ##outputs:
    #J is the matrix of partial derivatives (for chosen state and/or parameter derivative) 

    uvec=[u for _ in 1:nd+1] #creates a vector of vectors [x(t), x(t-τ_1),etc.]=[u,u,u,etc.]
    n=length(u) #length of vector x - the number of states (e.g. x1,x2,etc)
    m=length(f(uvec,pars)) #length of vector of outputs of f  e.g. if f is vector with 2 subfunctions, length of f =2

    ej=fill(0.0,n) #starts to set up unit vector

    #deepcopy used as function is handling vector of vectors
    u_plus=deepcopy(uvec) #ready for adding step to uvec (for state derivative)
    u_neg=deepcopy(uvec) #ready for subtracting step to uvec (for state derivative)
    pars_plus=deepcopy(pars) #ready for adding step to pars (for parameter derivative)
    pars_neg=deepcopy(pars) #ready for subtracting step to pars (for parameter derivative)

    if isempty(nx) && isempty(np)
        return println("Please state which derivative you would like")
    elseif length(nx)==1 && isempty(np) && isempty(v) #if given one state derivative to find 
        J=fill(NaN,m,n)
        for j in 1:n
            ej=fill(0.0,n)
            ej[j]=1
            u_plus[nx]=u_plus[nx].+h*ej #slight change (addition) to nx component of uvec [x(t), x(t-τ_1), etc.] -> e.g. if nx=1 you'd have [[u+h],[u],[u]]
            u_neg[nx]=u_neg[nx].-h*ej #slight change (subtraction) to nx component of the states

            J[:,j].=(1/2h)*(f(u_plus,pars) - f(u_neg,pars)) #finds finite difference with adjusted state values
            
            u_plus[nx]=deepcopy(uvec[nx]) #resets n_plus to what it was before (so it doesn't affect the next interation)
            u_neg[nx]=deepcopy(uvec[nx]) #resets n_neg to what it was before (so it doesn't affect the next interation)
        end 
        return J
    elseif isempty(nx) && length(np)==1 && isempty(v) #if given parameter derivative to find
        J=fill(NaN,n)
        pars_plus[np]=pars_plus[np]+h #slight change (addition) to the parameter at the index specified in np
        pars_neg[np]=pars_neg[np]-h ##slight change (subtraction) to the parameter at the index specified in np
        J=(1/2h)*(f(uvec,pars_plus)-f(uvec,pars_neg)) #finds finite difference with adjusted parameter values
            

        pars_plus[np]=deepcopy(pars[np]) #resets pars_plus to what it was before (so it doesn't affect the next interation)
        pars_neg[np]=deepcopy(pars[np]) #resets pars_plus to what it was before (so it doesn't affect the next interation)
        return J
    elseif length(nx)==1 && length(np)==1 && isempty(v) #if given state and parameter derivative to find
        J=fill(NaN,m,n)
        for j in 1:n
            ej=fill(0.0,n)
            ej[j]=1
            u_plus[nx]=u_plus[nx]+h*ej #slight change (by h) to state for derivative that was asked for 
            u_neg[nx]=u_neg[nx]-h*ej

            pars_plus[np]=pars_plus[np]+k #slight change (by k) to parameter indexed by np
            pars_neg[np]=pars_neg[np]-k

            J[:,j]=(1/(4*h*k))*(f(u_plus,pars_plus)-f(u_plus,pars_neg)-f(u_neg,pars_plus)+f(u_neg,pars_neg)) #finds finite difference with adjusted state and parameter values
            
            u_plus[nx]=deepcopy(uvec[nx]) #resets n_plus to what it was before (so it doesn't affect the next interation)
            u_neg[nx]=deepcopy(uvec[nx]) #resets n_neg to what it was before (so it doesn't affect the next interation)
            pars_plus[np]=deepcopy(pars[np])
            pars_neg[np]=deepcopy(pars[np])
        end
        return J
    end
end 