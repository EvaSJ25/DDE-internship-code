function j_eval(ti,te;diff=0)#j_eval(ti,te,n;diff=0)#, wj=[]) #(n) #n not needed as n=length(ti)
    ##inputs:
    #ti=the interpolation points/node
    #te=the points you want to evaluate at
    #diff= what derivative you want (0=evaluation, 1= 1st derivative,2= 2nd derivative )
    ##CURRENTLY FUNCTION ONLY WORKS FOR DIFF=0!!!!

    ##output:
    #

    nint=length(ti) #number of interpolation points/nodes (n+1)???? e.g. nint=3 => (x_0,x_1,x_2) so n=2
    nnew=length(te)#number of points wanted to be evaluated
    #lx=fill(NaN,nnew)#creates vector for l(x) values

    #if diff==0
    #elseif diff==1
    #elseif diff==2
    #end  

    #Computing the weights w_j
    #In the theory j runs from 0 to n for coding we'll run from 1 to n+1 (for Julia indices)
    wjvec=fill(NaN,nint)
    wjvec[1]=1 #w_0^(0)=1
    wjvec=fill(1.0,nint)
    for j in 2:nint #same as for j=1 to j=n in barycentric interpolation source
        for k in 1:j-1 #k in 0:j-1
            wjvec[k]=(ti[k]-ti[j])*wjvec[k]
        end 
        for t in ti[1:j-1]
            wjvec[j]*=(ti[j]-t)
        end 
    end 

    for j in 1:nint 
        wjvec[j]=1/wjvec[j]
    end 


    #l_j(x)=(wj/(x-x_j))/sum(wk/(x-x_k))

    numer=0.0
    denom=0.0
    lvals=fill(NaN,nnew,nint) #the values of l_j(x^(i))
    for i in 1:nnew
        for j in 1:nint
            numer=wjvec[j]/(te[i]-ti[j])
            #@infiltrate
            for k in 1:nint
                denom+=(wjvec[k]/(te[i]-ti[k]))
                #@infiltrate
            end 
            #lvals[i,j]=(lx[i]*wjvec[j])/(te[i]-ti[j])
            @infiltrate
            lvals[i,j]=numer/denom
            numer=0.0
            denom=0.0
        end 
    end 

    #return lx,wjvec,lvals
    #return lvals
    #return wjvec
    return lvals
end 