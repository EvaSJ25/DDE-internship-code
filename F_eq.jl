using LinearAlgebra
function F_eq(f_DDE, u0, pars,nd; par_indx=[]) #need to add delay equation later!!!
    #f_DDE is the RHS of the system
    #u0 is the initial guess (for the states x_i's)
    #pars are parameters of the system (including the delays)
    #nd is the number of delays
    #par_indx is the index of the parameter the user wishes to vary (to use continuation over)
    #Set par_indx=0 is you want to use the function output in Newton
    #p0 is the initial guess of the (varying) parameter; default is empty 
    
    n=length(u0) # number of states (x_i's)
    u0vec=[u0 for _ in 1:nd+1] #repeats the initial guess over all (delayed) states
    #nf=length(f_DDE(u0vec,pars)) #dimension of f
    np=length(par_indx)
    y0=vcat(u0,pars[par_indx]) #n+1 dimensional (n+1th term is the parameter being varied)


    #This works for use in newton and track_curve
    function f_equilibrium(y) #creates a function that can be used in newton and track_curve functions
    if par_indx==[]
        x=y[1:n]
        u=[x for _ in 1:nd+1]
    else  
        x,p=y[1:n],y[n+1]
        pars[par_indx]=p
        u=[x for _ in 1:nd+1]
    end 
    return f_DDE(u,pars) 
    end


    return y0,f_equilibrium
end 
