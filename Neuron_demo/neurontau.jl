function neurontau(uvec, pars) #delay functions for the 3 delays in neuron example
    #pars (the parameters) in this case are [κ, β, a12, a21,tau1, tau2, taus]
    tau1=pars[5] #5th given parameter is delay τ_1
    tau2=pars[6] #6th given parameter is delay τ_2
    taus=pars[7] #7th given parameter is delay τ_s
    tau=[tau1,tau2,taus] #vector of the 3 delays
    return tau
end 