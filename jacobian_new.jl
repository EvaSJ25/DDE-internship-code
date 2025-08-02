function jacobian_new(f,x; h=1e-5)
    n=length(x) #length of vector x
    #@infiltrate
    m=length(f(x)) #length of vector of outputs of f  e.g. if f is vector with 2 fubfuncton, length of f =2
    J=fill(NaN,m,n) #array for jacobian matrix J is a matrix with m rows and n columns 
    ej=fill(0.0,n) #defines ej vectors

    for j in 1:n
        ej=fill(0.0,n) #resets ej for each j iteration
        ej[j]=1 #sets jth component to 1 (all the others are zero)
        J[:,j]=(1/2h)*(f(x+(h*ej)) - f(x-(h*ej))) #differentiates function f with respect to x j-th component
    end

    return J
end