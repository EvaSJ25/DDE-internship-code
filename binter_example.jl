function binter_example(x;diff=0)
    if diff==0
        fun=abs.(x)+0.5*x -x.^2
        return fun
    elseif diff==1 
        fund1=x./abs.(x) -2*x .+ 0.5
        return fund1
    elseif diff==2
        fund2=-2
        return fund2
    end 
end