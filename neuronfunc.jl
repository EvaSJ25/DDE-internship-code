function neuronfunc(xx::Vector{Vector{Float64}},pars::Vector{Float64})
    #xx are the states (and delays) xx[1][2]=x2(t), x[2][1]=x1(t-tau1)
    #In this example pars=[kappa, betam a12, a21]
    x1dot=-pars[1]*xx[1][1] + pars[2]*tanh(xx[4][1]) + pars[3]*tanh(xx[3][2])
    x2dot=-pars[1]*xx[1][2] + pars[2]*tanh(xx[4][2]) + pars[4]*tanh(xx[2][1])  
    return [x1dot,x2dot]
end