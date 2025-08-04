using LinearAlgebra
function track_curve_new(userf,y0,ytan; h=1e-5,userdf=x->jacobian_new(userf,x;h=h), nmax=100,stepsize=0.01,tol=1e-8,maxit=6)
    n=length(y0) 
    m=length(ytan)
    s=stepsize
    ylist= [fill(NaN,n) for _ in 1:nmax+1] #creates array 
    ytan=ytan 
    ylist[1]=y0 #sets first component of ylist equal to inital guess
    en1=fill(0.0,n) #n basis vector (all components 0 except nth component)
    en1[n]=1 #assigns nth component of nth basis vector to 1

    for k in 2:nmax+1
        #@infiltrate
        ypred= ylist[k-1] +s.*ytan

        #@infiltrate
        #Below finds solution yk
        fbound = y-> vcat(userf(y), (ytan'*(y-ypred)))   
        #@infiltrate 
        newt=newton_new(fbound,ypred) #newton outputs
        ylist[k]=newt[1] #assigns solution yk to be equal to 1st output by newton
        converged=newt[2] #assigns converged to be second output of newton
        #@infiltrate
        if converged==false #if newton hasn't converged it will break the loop
            println("Has not converged") #tells user why the loop was broken (the solution didn't converge)
            break
        end 

        #Below finds new tangent
        A=vcat(userdf(ylist[k]),ytan')
        #@infiltrate
        z=A\en1 
        ytan=(z/norm(z,Inf))*(z'*ytan) #infinity norm used and new tangent found 
        #@infiltrate
    end
    return ylist, ytan #returns values along the curve for f and returns final tangent
end