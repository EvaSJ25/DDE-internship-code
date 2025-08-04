function tau3delays(uvec, pars) #3 delays for a 2 parameter system
    tau1=pars[3]
    tau2=pars[4]
    tau3=pars[5]
    tau=[tau1,tau2,tau3]
    return tau
end 