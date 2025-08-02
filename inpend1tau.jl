function inpend1tau(uvec,pars) #takes inputs: x(t), x(t-tau_1),..,x(t-tau_nd),parameters
    tau=pars[3] #for inverted pendulum example with one (constant) delay, tau is the 3rd parameter [a,b,tau]
    return tau
end 