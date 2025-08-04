#Trial with 2 delays
#xdot=v(t)
#vdot=sin(x(t))-ax(t-tau1)-bv(t-tau2)

function invpend_2delay(u,pars)
    xdot=u[1][2]
    vdot=sin(u[1][1]) - pars[1]*u[2][1]-pars[2]*u[3][2] #a=1, b=0.5 here
    return [xdot,vdot]
end