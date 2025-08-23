function j_diff(tbase)
    #input:
    #tbase are the interpolation points/nodes

    #output:
    # 1st-order differentiation matrix for the interpolation points

    nint=length(tbase)

    wjvec=fill(NaN,nint)
    wjvec[1]=1 #w_0^(0)=1
    wjvec=fill(1.0,nint) #set all weights to 1 to start
    for j in 2:nint #same as for j=1 to j=n in barycentric interpolation source
        for k in 1:j-1 #k in 0:j-1
            wjvec[k]=(tbase[k]-tbase[j])*wjvec[k]
        end 
        for t in tbase[1:j-1]
            wjvec[j]*=(tbase[j]-t)
        end 
    end 

    for j in 1:nint
        wjvec[j]=1/wjvec[j] 
    end 

    D1=fill(0.0,nint,nint) #creates blank array for 1st-order differentiation matrix
    for i in 1:nint
        for j in 1:nint
            if i!=j
                D1[i,j]=(wjvec[j])/(wjvec[i]*(tbase[i]-tbase[j])) 
            end 
        end
        D1[i,i]=-sum(D1[i,:])
    end 

    return D1
end 