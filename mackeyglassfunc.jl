function mackeyglassfunc(u,pars)
    #pars=[beta,gamma,n] (e.g.=[1.0,1,2]) (pars components need to be in Float64 form)
    #u=[[x(t)],[x(t-tau)]] (components need to in vector (and Float64) form)
    #@infiltrate
    u1=convert(Float64,u[1][1])
    u2=convert(Float64,u[2][1])
    #xdot=pars[1]*((u2)/(1+u2^Int(pars[3])))-pars[2]*u1
    #xdot=pars[1]*((u2)/(1+real((u2+0im)^(pars[3]))))-pars[2]*u1
    xdot=pars[1]*((u2)/1+(real((u2+0im)^(pars[3]))))-pars[2]*u1

    #@infiltrate

    return [xdot]
end 