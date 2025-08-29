using LinearAlgebra
function create_ststfunc(f_DDE, u0, pars,nd; par_indx=[])
    ##inputs:
    #f_DDE is the RHS of the system
    #u0 is the initial guess (for the states x_i's)
    #pars are parameters of the system (including the delays), ensure that p0 - initial guess of the (varying) parameter - is in this vector
    #nd is the number of delays
    #par_indx is the index of the parameter the user wishes to vary (to use continuation over)
    #Leave par_indx empty to use the function output in Newton function

    ##outputs:
    #y0 is a vector containing the initial guess u0 and initial guess of parameter that u0 is an equilibrium/steady-state for
    #f_equilibrium is a function that is used to find teh equilibrium (finds f([u,u,..,u],pars))

    n=length(u0) # number of states (x_i's)
    u0vec=[u0 for _ in 1:nd+1] #repeats the initial guess over all (delayed) states
    np=length(par_indx) #number of varied parameters (e.g. 0 for use in newton function, 1 for finding equilibrium)
    y0=vcat(u0,pars[par_indx]) #n+1 dimensional vector of initial guess for equilibirum information (u0 with the n+1th term being the parameter being varied))

    #This works for use in newton and track_curve
    function f_equilibrium(y) #creates a function that can be used in newton and track_curve functions to find equilibria
    if par_indx==[] #for use in newton function
        x=y[1:n]
        u=[x for _ in 1:nd+1]
    else  #for use in track_curve function
        x,p=y[1:n],y[n+1]
        pars[par_indx]=p
        u=[x for _ in 1:nd+1]
    end 
    return f_DDE(u,pars) 
    end

    return y0,f_equilibrium
end 
