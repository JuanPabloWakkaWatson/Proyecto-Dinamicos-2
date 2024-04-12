using DifferentialEquations
using Plots
using LaTeXStrings

# Definir el sistema de EDOs para los osciladores acoplados
function osciladores!(du, u, p, t)
    x1, x2 = u
    t1, t2, I1, I2 = p
    du[1] = -x1/t1 + I1  # Ecuación para x1
    du[2] = -x2/t2 + I2  # Ecuación para x2
end

# Parámetros del sistema: t1, t2, I1, I2
p = [3.0, 2.0, 2.0 , -3.0]  # Valores para t1, t2, I1, I2

# Condiciones iniciales
u0 = [0.5, 0.6]

# Rango de tiempo
tspan = (0.0, 50.0)

# Resolver el sistema de EDOs
prob = ODEProblem(osciladores!, u0, tspan, p)
sol = solve(prob)

# Graficar solución integral (x1 y x2 vs t)
p1 = plot(sol, vars=(1), xlabel=L"t", ylabel=L"x_1", label="x_1", title="Solución Integral vs Tiempo")
p2 = plot(sol, vars=(2), xlabel=L"t", ylabel=L"x_2", label="x_2", title="Solución Integral vs Tiempo")
plot(p1, p2, layout=(2, 1), legend=:topright)

# Graficar x1 vs x2
plot(sol[1,:], sol[2,:], title="Campo de fases", xlabel=L"x_1", ylabel=L"x_2", label="", legend=false)
scatter!([u0[1]], [u0[2]], color=:red, markersize=5, markerstrokewidth=0, label="")
