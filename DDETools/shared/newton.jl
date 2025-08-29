using LinearAlgebra
function newton(f,x0; h=1e-5,df=x->jacobian(f,x;h=h),tol=1e-8,maxit=100)
    ##inputs:
    #f is the RHS of system 
    #x0 is starting guess for x
    #h is the stepsize needed for jacobian function 
    #tol is the tolerance for which the two convergence conditions must be fulfilled
    #maxit is the number of maximum iterations the function is allowed to do (if the function does not fulfil the tolerance conditions and exceeds maxit then the function did not converge succesfully onto a point x)

    ##outputs:
    #x is the value of x that fulfils newton (if converged)
    #converged tells user if their run was successful and x has converged to a point
    #J is the Jacobian at the previous point to the outputted x

    n=length(x0) #dimensions of x
    x=x0-df(x0)\f(x0) #calculates x for first step

    xold=x0 #makes first xold (previous step) equal to x0
    iter=0 #iter is the iteration step
    converged=(norm(x-xold) < tol) && (norm(f(x))<tol) #conditions for convergence
    
    #While both the convergence conditions are not fulfilled and while the iteration number is below the maximum iteration number the cycle continues
    #Note: once the maximum iteration is reached the final value of x found is returned
    while (norm(x-xold) > tol || norm(f(xold))>tol) && (iter<maxit) 
        xold=x #reinitialises xold as the previously calculated x
        x=xold - (df(xold)\f(xold)) #finds new point x
        converged = (norm(x-xold) < tol) && (norm(f(x))<tol) #converged="true" if both conditions with the tolerance are met
        iter+=1 #index increased by 1 while the loop is still continued (the conditions are still being fulfilled, converged ="false")
    end
    J=df(xold) #jacobian at previous x 
    return x,converged, J
end