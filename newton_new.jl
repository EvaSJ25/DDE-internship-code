using LinearAlgebra
function newton_new(f,x0; h=1e-5,df=x->jacobian_new(f,x;h=h),tol=1e-8,maxit=100)#maxit increased from 6 to 100
    n=length(x0)
    x=x0-df(x0)\f(x0) #calculates x for first step

    xold=x0 #makes first xold (previous step) equal to x0
    iter=0 #iter is the iteration step
    converged=(norm(x-xold) < tol) && (norm(f(x))<tol)
    
    #While the difference between the previous and new values of x and f(xold) are greater than some tolerance the cycle continues
    #Also while the number of iterations is less than a set maximum (once the maximum is reached the final value of x found is returnesd)
    while (norm(x-xold) > tol || norm(f(xold))>tol) && (iter<maxit) 
        xold=x
        x=xold - (df(xold)\f(xold))
        converged = (norm(x-xold) < tol) && (norm(f(x))<tol) #converged="true" if both conditions with the tolerance are met
        iter+=1 #index increased by 1 while the loop is still continued (the conditions are still being fulfilled)
    end
    J=df(xold) #jacobian at previous x 
    return x,converged, J
end