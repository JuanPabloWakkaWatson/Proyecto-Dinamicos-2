using DifferentialEquations
using Plots
using LaTeXStrings

include("fieldlorenz.jl")

x0 = 0.6
y0 = 0.5
z0 = 0.1
u0 = [x0,y0,z0]
p = (10,28,8/3)
tspan = (0.0,50.0)

prob = ODEProblem(fieldlorenz,u0,tspan,p)
sol = solve(prob)

plot(sol,vars=(1,2,3),xlabel=L"x",ylabel=L"y",zlabel=L"z",label="")
scatter!([x0],[y0],[z0],color=:red,markersize=5,markerstrokewidth=0,label="")