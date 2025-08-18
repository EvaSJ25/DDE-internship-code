function j_eval_multidim(ti,te,fvec;diff=0)#(f,ti,te)#used for when x is multidimensional (used for an dimension of x??)
    #doesn't work for 1D x yet (problem with fvec[i,j] in line 78ish - because j not existent for 1D x)
    #f is the (known) function - only useful if you know it
    #ti=the interpolation points/node
    #te=the points you want to evaluate at
    #fvec are the values of f(x1j,x2i,...,xndim i) #f_ij=f(xi,yj) #probably in matrix form!
    #diff= what derivative you want (0=evaluation, 1= 1st derivative,2= 2nd derivative )

    ##ouput:
    #lvals=  #J=jacobian
    
    
    #for val in te if val in ti
        #fx[]=fvec[]
    #else 
        #...
    #end
    #end 
    ndim=map(length,ti)[1]#finds the dimension of x (number of states)
    #Currently we're assuming nx=ny=n 
    nint=length(ti) #number of interpolation points/nodes (n+1)???? e.g. nint=3 => (x_0,x_1,x_2) so n=2
    nnew=length(te)#number of points wanted to be evaluated
    #nnew=map(length,te)[1]#number of points wanted to be evaluated

    lx=fill(NaN,nnew)#creates matrix for l(x) values
    tivals=[fill(NaN,nint) for _ in 1:ndim]
    tevals=[fill(NaN,nnew) for _ in 1:ndim]
    #@infiltrate
    for k in 1:ndim
        tivals[k]=[u[k] for u in ti]#ti_k values #interpoltaion points t1 values for k=1, t2 values for k=2, etc.
        #@infiltrate
        tevals[k]=[u[k] for u in te]#kth row corresponds to.. (e.g for ndim=2 (x,y), k=1 corresponds to x, k=2 corresponds to y values)
    end 

    #@infiltrate
    #for i in 1:nnew
     #   lx[i]=1
      #  for t_i in ti
       #     lx[i]*=(te[i]-t_i)
            #@infiltrate
        #end
   # end 

    #Computing the weights w_j
    #In the theory j runs from 0 to n for coding we'll run from 1 to n+1 (for Julia indices)
    wjvec=fill(NaN,nint,ndim) #matrix where the rows correspond to j and the columns to dimension of x (e.g. col 1 is x, col 2 is y) - i.e wjvec[1,1]=x_j^(x)
    #wjvec[1]=1 #w_0^(0)=1
    wjvec=fill(1.0,nint,ndim) #set all weights to begin at 1
    for i in 1:ndim
        for j in 2:nint #same as for j=1 to j=n in barycentric interpolation source
            for k in 1:j-1 #k in 0:j-1
                #wjvec[i,k]=(tivals[i][k]-tivals[i][j])*wjvec[i,k]
                wjvec[k,i]=(tivals[i][k]-tivals[i][j])*wjvec[k,i]

            end 
            for t in tivals[i][1:j-1]
                #@infiltrate
                #wjvec[i,j]*=(tivals[i][j]-t)
                wjvec[j,i]*=(tivals[i][j]-t)

            end 
        end 
        for j in 1:nint 
        #wjvec[i][j]=1/wjvec[i][j]
        wjvec[j,i]=1/wjvec[j,i]

    end
    end 

    #wjvec columns correspond to dim of x so if you want x out of (x,y) look at column 1! (these are your w_i's)

    gvec=fill(NaN,nint,nnew) #creates (matrix) array for g_i(y) = col k is for new points (x_k,y_k)
    numer=0.0 #set starting numerator as 0 
    denom=0.0 #set starting denominator as 0
    for k in 1:nnew
        for i in 1:nint
            for j in 1:nint
                numer+=((wjvec[j,2])/(tevals[2][k]-tivals[2][j]))*fvec[i,j]#fvec[i,k] #wjvec[,2] is for y need to sort it out for 3+ dimensional x!!!!
                denom+=((wjvec[j,2])/(tevals[2][k]-tivals[2][j]))
                #@infiltrate
            end 
            gvec[i,k]=numer/denom
            #@infiltrate
            numer=0.0 #reset numerator
            denom=0.0 #reset denominator
        end 
    end 


    pvals=fill(NaN,nnew)#,ndim) can have fill(NaN,nnew,nf) where nf is the dimension of f(x)
    numer2=0.0
    denom2=0.0
    for k in 1:nnew
        for i in 1:nint
            numer2+=((wjvec[i,1])/(tevals[1][k]-tivals[1][i]))*gvec[i,k]#fvec[i,k] #wjvec[,1] is for x - need to sort it out for 3+ dimensional x!!!!
            denom2+=((wjvec[i,1])/(tevals[1][k]-tivals[1][i]))
            
            pvals[k]=numer2/denom2
            #pvals[k,i]=numer2/denom2
            #@infiltrate
        end 
        pvals[k]
        numer2=0.0
        denom2=0.0
    end 


   # lvals=fill(NaN,nnew,nint) #the values of l_j(x^(i))
    #for i in 1:nnew
     #   for j in 1:nint
      #      lvals[i,j]=(lx[i]*wjvec[j])/(te[i]-ti[j])
       # end 
    #end 

    #return lx,wjvec,lvals
    #return lvals
    #return wjvec
    return pvals
end 