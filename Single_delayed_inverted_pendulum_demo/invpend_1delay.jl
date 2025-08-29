#Inverted pendulum with 1 delay, tau
#xdot=v(t)
#vdot=sin(x(t))-ax(t-tau)-bv(t-tau)
function invpend_1delay(u::Vector{Vector{Float64}},pars::Vector{Float64})
    ##inputs:
    #u are the values of the states (x(t), x(t-τ)) in vector of vector form
    #pars are the system parameters

    #output: the values of the system for the given u and parameters
    xdot=u[1][2]
    vdot=sin(u[1][1]) - pars[1]*u[2][1]-pars[2]*u[2][2] #a=pars[1], b=pars[2]
    return [xdot,vdot]
end 