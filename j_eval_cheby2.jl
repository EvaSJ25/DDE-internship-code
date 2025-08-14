function j_eval_cheby2(ti,te;diff=0)#j_eval(ti,te,n;diff=0)#, wj=[]) #(n) #n not needed as n=length(ti)
    ##inputs:
    #ti=the interpolation points/node
    #te=the points you want to evaluate at
    #diff= what derivative you want (0=evaluation, 1= 1st derivative,2= 2nd derivative )

    ##ouput:
    #J=jacobian 

    nint=length(ti) #number of interpolation points/nodes (n+1)???? e.g. nint=3 => (x_0,x_1,x_2) so n=2
    nnew=length(te)#number of points wanted to be evaluated
    lx=fill(NaN,nnew)#creates vector for l(x) values

    for i in 1:nnew
        lx[i]=1
        for t_i in ti
            lx[i]*=(te[i]-t_i)
            #@infiltrate
        end
    end 

    #Computing the weights w_j
    #do i need to find w_j or should i just use Chebyshev points of the second kind?
    #In the theory j runs from 0 to n for coding we'll run from 1 to n+1 (to not confuse indices)
    #if wj==[]
        #w_0^(0)=1
     #   for j in 1:n 
      #      for k in 0:j-1
       #     end 
        #end 
        #for j in 0:n 
            #wj[j]=
        #end 
    #else
     #   wj=wj
    #end

    #To try and get a function working we'll take wj in form of Chebyshev points of the second kind:
    #wjvec=fill(NaN,n+1)
    wjvec=fill(NaN,nint)
    for j in 1:nint#j in 1:n+1
        if j==1 || j==nint#n+1  #if j=0 or j=n (in source)
            wjvec[j]=(-1)^(j-1)*(1/2)
        else 
            wjvec[j]=(-1)^(j-1)*1
        end
    end

    #@infiltrate
    lvals=fill(NaN,nnew,nint) #the values of l_j(x^(i))
    for i in 1:nnew
        for j in 1:nint
            lvals[i,j]=(lx[i]*wjvec[j])/(te[i]-ti[j])
            #@infiltrate
        end 
    end 

    #return lx,wjvec,lvals
    return lvals
end 