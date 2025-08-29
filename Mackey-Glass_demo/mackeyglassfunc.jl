#The Mackey-Glass equation:
#xdot=β[x(t-τ)/(1+x(t-τ)^n)]-γx(t)
function mackeyglassfunc(u,pars)
    ##inputs:
    #u=[[x(t)],[x(t-tau)]] (components need to be in Vector{Float64} form)
    #pars are the parameters = [beta,gamma,n] in Vector{Float64} form

    u1=convert(Float64,u[1][1]) 
    u2=convert(Float64,u[2][1])
    
    xdot=pars[1]*((u2)/(1+(real((u2+0im)^(pars[3])))))-pars[2]*u1

    return [xdot]
end 