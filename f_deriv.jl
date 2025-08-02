function f_deriv(f,u,pars,nd;nx=[],np=[],v=[],h=1e-5,k=1e-5) #u many columns
    #f is the system (equations)
    #u is the point jacobian being taken at
    #pars are the system parameters
    #nd is number of delays
    #nx is the state derivative you want (empty=no state derivative wanted)
    #nx=1 is derivtive wrt to state at time t, nx=2 is derivative wrt to state at t-tau_1
    #np is whether you want parameter derivative (empty = no parameter derivative wanted)
    #h is the small step for finite difference
    #k is the small steps for paramater (used only is nx and np both non-empty)
    uvec=[u for _ in 1:nd+1]
    n=length(u) #length of vector x - the number of states (e.g. x1,x2,etc)
    #@infiltrate
    m=length(f(uvec,pars)) #length of vector of outputs of f  e.g. if f is vector with 2 fubfuncton, length of f =2
    num_p=length(pars)

    ej=fill(0.0,n) 

    u_plus=deepcopy(uvec)
    u_neg=deepcopy(uvec)
    pars_plus=deepcopy(pars)
    pars_neg=deepcopy(pars)

    #if length(nx)==0 && length(np)==0
    if isempty(nx) && isempty(np)
        return println("Please state which derivative you would like")
    elseif length(nx)==1 && isempty(np) && isempty(v)
        J=fill(NaN,m,n)
        for j in 1:n
            ej=fill(0.0,n)
            ej[j]=1
            #@infiltrate
            u_plus[nx]=u_plus[nx].+h*ej
            #u_plus[nx]=u_plus[nx]+h*ej

            #u_plus[nx[1]]=u_plus[nx[1]]+h*ej #MATLAB version manual has this!

            u_neg[nx]=u_neg[nx].-h*ej
            #u_neg[nx]=u_neg[nx]-h*ej

            #@infiltrate
            J[:,j].=(1/2h)*(f(u_plus,pars) - f(u_neg,pars))
            #Below may not be needed
            u_plus[nx]=deepcopy(uvec[nx]) #resets n_plus to what it was before (so it doesn't affect the next interation)
            u_neg[nx]=deepcopy(uvec[nx])
        end 
        return J
    elseif isempty(nx) && length(np)==1 && isempty(v)
        J=fill(NaN,n)
        pars_plus[np]=pars_plus[np]+h
        J=(1/h)*(f(uvec,pars_plus)-f(uvec,pars))
        return J
    elseif length(nx)==1 && length(np)==1 && isempty(v)
        J=fill(NaN,m,n)
        for j in 1:n
            ej=fill(0.0,n)
            ej[j]=1
            u_plus[nx]=u_plus[nx]+h*ej
            u_neg[nx]=u_neg[nx]-h*ej

            pars_plus[np]=pars_plus[np]+k
            pars_neg[np]=pars_neg[np]-h
            #@infiltrate

            J[:,j]=(1/(4*h*k))*(f(u_plus,pars_plus)-f(u_plus,pars_neg)-f(u_neg,pars_plus)+f(u_neg,pars_neg))
            #@infiltrate
            u_plus[nx]=deepcopy(uvec[nx]) #resets n_plus to what it was before (so it doesn't affect the next interation)
            u_neg[nx]=deepcopy(uvec[nx])
            pars_plus[np]=deepcopy(pars[np])
            pars_neg[np]=deepcopy(pars[np])
        end
        return J
    elseif length(nx)==2 #&& ~isempty(v)
        J=fill(NaN,m,n)
        A=fill(NaN,m,n)#df/dx_nx[1]
        A1=fill(NaN,m,n)#df/dx_nx[2]        
        for k in 1:n 
            #for i in 1:2
            #u_plus[nx[i]]=u_plus[nx[i]]+h*ej
            #@infiltrate
            u_plus[nx[1]]=u_plus[nx[1]]+h*ej
            u_neg[nx[1]]=u_neg[nx[1]]-h*ej

            A[:,k]=(1/2h)*(f(u_plus,pars) - f(u_neg,pars))
            u_plus[nx[1]]=deepcopy(uvec[nx[1]]) #resets n_plus to what it was before (so it doesn't affect the next interation)
            u_neg[nx[1]]=deepcopy(uvec[nx[1]])

            #A[:,k]=(1/(h^2))*(f(u_plus,pars)-2*f(uvec,pars) + f(u_neg,pars))
            u_plus[nx[2]]=u_plus[nx[2]]+h*ej
            u_neg[nx[2]]=u_neg[nx[2]]-h*ej

            A1[:,k]=(1/2h)*(f(u_plus,pars) - f(u_neg,pars))
            u_plus[nx[2]]=deepcopy(uvec[nx[2]]) #resets n_plus to what it was before (so it doesn't affect the next interation)
            u_neg[nx[2]]=deepcopy(uvec[nx[2]])

            J=A1*(A*v)
        end 
        return J
    end
end 