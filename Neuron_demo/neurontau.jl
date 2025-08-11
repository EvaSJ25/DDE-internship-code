function neurontau(uvec, pars) #delay functions for the 3 delays in neuron example
    #pars (the parameters) in this case are [κ, β, a12, a21,tau1, tau2, taus]
    tau1=pars[5]
    tau2=pars[6]
    taus=pars[7]
    tau=[tau1,tau2,taus] #vector of the 3 delays
    return tau
end 