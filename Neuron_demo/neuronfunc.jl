#The functions of the system for the neuron example:
#x1dot=-κx1(t) + βtanh(x1(t-τ_s))+a12tanh(x2(t-τ_2))
#x2dot=-κx2(t) + βtanh(x2(t-τ_s))+a21tanh(x1(t-τ_1))
function neuronfunc(xx::Vector{Vector{Float64}},pars::Vector{Float64})
    ##inputs:
    #xx are the states (and delays) xx[1][2]=x_2(t), x[2][1]=x_1(t-tau1)
    #pars are the parameters =[κ, β, a12, a21,tau1, tau2, taus], where tau1, tau2, taus are the delays
    
    x1dot=-pars[1]*xx[1][1] + pars[2]*tanh(xx[4][1]) + pars[3]*tanh(xx[3][2])
    x2dot=-pars[1]*xx[1][2] + pars[2]*tanh(xx[4][2]) + pars[4]*tanh(xx[2][1])  
    return [x1dot,x2dot]
end