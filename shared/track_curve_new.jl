using LinearAlgebra
function track_curve_new(userf,y0,ytan; h=1e-5,userdf=x->jacobian_new(userf,x;h=h), nmax=100,stepsize=0.01,tol=1e-8,maxit=6)
    ##inputs:
    #userf is the function the user wishes to track (in a form such that it can be used in track curve - y0 is one dimension more than original x0 so the tangent can be used)
    #y0 is the initial point - it contains x0 (initial states) and a starting point for the varying parameter
    #ytan is the initial tangent
    #h is the stepsize for use in jacobian_new
    #tol is the tolerance used in newton_new
    #maxit are the maximum number of iterations newton_new is allowed to do
    #userdf is the jacobian_new function adapted for use in track_curve_new
    #nmax is the number of points you want to track (tracking more points gives a more complete visual of the tracked function)
    #stepsize is the stepsize added to the tangent (ytan) when predicting the new tracked point

    ##outputs:
    #ylist are all the points that have been tracked
    #ytan is the final tangent
    
    n=length(y0) 
    m=length(ytan)
    s=stepsize
    ylist= [fill(NaN,n) for _ in 1:nmax+1] #creates array 
    ytan=ytan 
    ylist[1]=y0 #sets first component of ylist equal to inital guess
    en1=fill(0.0,n) #n basis vector (all components 0 except nth component)
    en1[n]=1 #assigns nth component of nth basis vector to 1

    for k in 2:nmax+1
        ypred= ylist[k-1] +s.*ytan #finds prediction for next ylist entry

        #Below finds solution yk
        fbound = y-> vcat(userf(y), (ytan'*(y-ypred)))   
        newt=newton_new(fbound,ypred) #newton outputs
        ylist[k]=newt[1] #assigns solution yk to be equal to 1st output by newton
        converged=newt[2] #assigns converged to be second output of newton
        
        if converged==false #if newton hasn't converged it will break the loop
            println("Has not converged") #tells user why the loop was broken (the solution didn't converge)
            break
        end 

        #Below finds new tangent
        A=vcat(userdf(ylist[k]),ytan')
        z=A\en1 
        ytan=(z/norm(z,Inf))*(z'*ytan) #infinity norm used and new (normalised) tangent found 
    end
    return ylist, ytan #returns values along the curve for f and returns final tangent
end