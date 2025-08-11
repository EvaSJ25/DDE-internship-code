function inpend2tau(uvec,pars)
    ##inputs: 
    #uvec is the vector of vectors for x(t), x(t-tau_1),..,x(t-tau_nd)
    #pars are the parameters

    #outputs: tau for 2 delayed inverted pendulum example 
    tau=[pars[3],pars[4]]
    return tau
end 