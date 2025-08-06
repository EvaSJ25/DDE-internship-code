#Trial with 2 delays
#xdot=v(t)
#vdot=sin(x(t))-ax(t-tau1)-bv(t-tau2)

function invpend_2delay(u,pars)
    ##inputs:
    #u are the values of the states (x(t), x(t-τ)) in vector of vector form
    #pars are the system parameters

    #output: the values of the system for the given u and parameters
    xdot=u[1][2]
    vdot=sin(u[1][1]) - pars[1]*u[2][1]-pars[2]*u[3][2] #a=1, b=0.5 here
    return [xdot,vdot]
end