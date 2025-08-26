using LinearAlgebra
function newton(f,x0; h=1e-5,df=x->jacobian(f,x;h=h),tol=1e-8,maxit=100)
    ##inputs:
    #f is the RHS of system 
    #x0 is starting guess for x
    #h is the stepsize needed for jacobian_new
    #tol is the tolerance for which the two convergence conditions must be fulfilled
    #maxit is the number of maximum iterations the function is allowed to do (if the function does not fulfil the tolerance conditions and exceeds maxit thenthe function did not converge succesfully onto a point x)

    ##outputs:
    #value of x that fulfils newton (if converged)
    #converged tells user if their run was successful and x has converged to a point
    #J is the Jacobian at the previous point of the outputted x
    n=length(x0)
    x=x0-df(x0)\f(x0) #calculates x for first step

    xold=x0 #makes first xold (previous step) equal to x0
    iter=0 #iter is the iteration step
    converged=(norm(x-xold) < tol) && (norm(f(x))<tol) #conditions for convergence
    
    #While the difference between the previous and new values of x and f(xold) are greater than some tolerance the cycle continues
    #Also while the number of iterations is less than a set maximum (once the maximum is reached the final value of x found is returnesd)
    while (norm(x-xold) > tol || norm(f(xold))>tol) && (iter<maxit) 
        xold=x #reinitialises xold as the perviously calculated x
        x=xold - (df(xold)\f(xold)) #finds new point x
        converged = (norm(x-xold) < tol) && (norm(f(x))<tol) #converged="true" if both conditions with the tolerance are met
        iter+=1 #index increased by 1 while the loop is still continued (the conditions are still being fulfilled)
    end
    J=df(xold) #jacobian at previous x 
    return x,converged, J
end