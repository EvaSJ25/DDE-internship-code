function inpend1tau(uvec,pars) 
    ##inputs: 
    #uvec is the vector of vectors for x(t), x(t-tau_1),..,x(t-tau_nd)
    #pars are the parameters

    #outputs: tau
    
    tau=pars[3] #for inverted pendulum example with one (constant) delay, tau is the 3rd parameter [a,b,tau]
    return tau
end 