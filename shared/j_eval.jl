function j_eval(ti,te;diff=0)#j_eval(ti,te,n;diff=0)#, wj=[]) #(n) #n not needed as n=length(ti)
    ##inputs:
    #ti=the interpolation points/node
    #te=the points you want to evaluate at
    #diff= what derivative you want to find points (te) at (0=evaluation, 1= 1st derivative,2= 2nd derivative )

    ##output:
    #Dk is matrix such that when multiplied by fvec (the values of f(ti)), the user gets the values of f (when diff=0), f' (when diff=1), etc. at the desired points, te

    nint=length(ti) #number of interpolation points/nodes
    nnew=length(te)#number of points wanted to be evaluated

    #Computing the weights w_j
    #In the theory (Berrut et al  2004), j runs from 0 to n for coding we'll run from 1 to n+1 (for Julia indices)
    wjvec=fill(NaN,nint)
    wjvec[1]=1 #w_0^(0)=1
    wjvec=fill(1.0,nint) #start all weights equal to 1
    for j in 2:nint 
        for k in 1:j-1 
            wjvec[k]=(ti[k]-ti[j])*wjvec[k]
        end 
        for t in ti[1:j-1]
            wjvec[j]*=(ti[j]-t)
        end 
    end 

    for j in 1:nint
        wjvec[j]=1/wjvec[j] 
    end 

    numer=0.0
    denom=0.0
    lvals=fill(NaN,nnew,nint) #the values of l_j(t^(i))
    for i in 1:nnew
        for j in 1:nint
            if te[i]==ti[j]
                lvals[i,j]=1.0 # to avoid dividing by Inf (for xi=xj, then p(xi)=f(xj))
            else
                numer=wjvec[j]/(te[i]-ti[j])
                for k in 1:nint
                    denom+=(wjvec[k]/(te[i]-ti[k]))
                end 
                lvals[i,j]=numer/denom #This formula is obtained from Equation (4.2) in (Berrut et al 2004)
                numer=0.0 #resets numerator for next iteration
                denom=0.0 #resets denominator for next iteration
            end
        end 
    end 

    E=lvals #matrix of l_j(t_i) values
    D=j_diff(ti) #Finds 1st-order derivative for interpolation points (D defined as D1 in (Berrut et al 2004))

    Dk=fill(NaN,nnew, nint) #blank arry for matrix to multiply with fj vector to get derivatives of te for interpolated function
    d=diff #number of derivatives to be taken

    for d=diff
        Dk=E*D^(d)
    end 
    return Dk

    #Reference: Barycentric Lagrange Interpolation, J-P. Berrut, L.N. Trefethen, SIAM Review, Vol 46, No.3 pp. 501-517, 2004
end 