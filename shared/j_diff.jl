function j_diff(tbase)
    #input:
    #tbase are the interpolation points/nodes

    #output:
    #D1 is the 1st-order differentiation matrix for the interpolation points (D_{ij}= l'j(t_i) for an interpolation point t_i

    nint=length(tbase) #number of interpolation points

    wjvec=fill(NaN,nint) #creates blank array for weights (wj)
    wjvec[1]=1 #w_0^(0)=1
    wjvec=fill(1.0,nint) #set all weights to 1 to start

    for j in 2:nint 
        for k in 1:j-1 
            wjvec[k]=(tbase[k]-tbase[j])*wjvec[k] #updates weight with next interpolation point
        end 
        for t in tbase[1:j-1]
            wjvec[j]*=(tbase[j]-t) #updates weight j 
        end 
    end 

    for j in 1:nint
        wjvec[j]=1/wjvec[j] #divides not to avoid extra unnecessary divisions (wj formula given in Equation (3.2) in (Berrut et al 2004))
    end 

    D1=fill(0.0,nint,nint) #creates blank array for 1st-order differentiation matrix
    for i in 1:nint
        for j in 1:nint
            if i!=j
                D1[i,j]=(wjvec[j])/(wjvec[i]*(tbase[i]-tbase[j])) #Equation (9.4) in (Berrut et al 2004)
            end 
        end
        D1[i,i]=-sum(D1[i,:]) #Equation (9.5) in (Berrut et al 2004)
    end 
    return D1
    #Note: Formula in above is taken from -> Barycentric Lagrange Interpolation, J-P. Berrut, L.N. Trefethen, SIAM Review, Vol 46, No.3 pp. 501-517, 2004 
end 