function fieldlorenz(du,u,p,t)
    a,b,c = p
    du[1] = a*(u[2]-u[1])
    du[2] = u[1]*(b-u[3]) - u[2]
    du[3] = u[1]*u[2] - c*u[3]
end