using DifferentialEquations
using Plots
using LaTeXStrings

# Función gaussiana para modelar las interacciones
function gaussian(t, mean, std)
    1 / (std * sqrt(2 * π)) * exp(-0.5 * ((t - mean) / std)^2)
end

# Usar destructuración en la función osciladores!
function osciladores!(du, u, p, t)
    x1, x2 = u
    τ1, τ2, I1, I2, ε, firing_times1, firing_times2, std_g = p
    
    # Calcular E1 y E2 usando gaussianas
    E1 = sum(gaussian(t, tf + 1.0, std_g) for tf in firing_times1)  # Retardo de 1.0 para E1
    E2 = sum(gaussian(t, tf, std_g) for tf in firing_times2)

    # Ecuaciones con términos de acoplamiento
    du[1] = -x1/τ1 + I1 + ε * E2  # Agregar E2 a la dinámica de x1
    du[2] = -x2/τ2 + I2 + ε * E1  # Agregar E1 a la dinámica de x2
end

params = (
    τ1 = 3.0,
    τ2 = 2.0,
    I1 = 2.0,
    I2 = -3.0,
    ε = 0.5,
    firing_times1 = [5, 15, 25, 35, 45],
    firing_times2 = [10, 20, 30, 40, 50],
    std_g = 0.1
)

# Condiciones iniciales y rango de tiempo
u0 = [0.5, 0.6]
tspan = (0.0, 50.0)

# Resolver el sistema de EDOs
prob = ODEProblem(osciladores!, u0, tspan, params)
sol = solve(prob, reltol=1e-6, abstol=1e-6)

# Graficar solución integral
plot(sol, xlabel="t", label=["x1" "x2"], layout=(2, 1), legend=:topright)

# Graficar el campo de fases
p_phase = plot(title="Retrato fase", xlabel=L"x_1", ylabel=L"x_2", legend=false)
plot!(p_phase, sol[1,:], sol[2,:], label="", color=:blue)
scatter!(p_phase, [u0[1]], [u0[2]], label="Inicio", color=:red, markersize=5)

# Mostrar todas las gráficas
plot(p_phase)
