module DDETools

include("shared/jacobian.jl") #Jacobian function 
include("shared/newton.jl") #Newton function
include("shared/track_curve.jl") #Tracking curve function (used for continuation and for finding equilibria branches)
include("shared/f_deriv.jl") #used to find state and/or parameter derivatives
include("shared/create_ststfunc.jl") #Initialisates a function to find equilibria 
include("shared/stab_func_matrix.jl") #uses large matrix method to approximate the stability and find the eigenvalues (and eigenvectors and ω if user interested in Hopf bifurcation)
include("shared/j_diff.jl") #finds the differentiation matrix of all interpolation points x_j
include("shared/j_eval.jl") #uses barycentric interpolation to find the l_j(xi) (Lagrange polynomial) values needed for use in stab_func_DDE
include("shared/stab_func_DDE.jl") #uses method explored in (Breda et al 2009) to find stability and eigenvalues (and if asked, also eigenvectors) of equilibria
include("shared/create_hopffunc.jl") #finds initial guess for Hopf information and creates a function that can find the true Hopf values
include("shared/create_foldfunc.jl") #finds initial guess for fold information and creates a function that can find the true fold values

end