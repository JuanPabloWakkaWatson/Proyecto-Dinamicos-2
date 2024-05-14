using DifferentialEquations
using Plots
using LaTeXStrings

# Función gaussiana para modelar las interacciones
function gaussian(t, mean, std)
    1 / (std * sqrt(2 * π)) * exp(-0.5 * ((t - mean) / std)^2)
end

# Función para las ecuaciones diferenciales
function osciladores!(du, u, p, t)
    x1, x2 = u
    τ1, τ2, I1, I2, ε, firing_times1, firing_times2, std_g = p
    
    # Calcular E1 y E2 usando gaussianas
    E1 = sum(gaussian(t, tf + 1.0, std_g) for tf in firing_times1)
    E2 = sum(gaussian(t, tf, std_g) for tf in firing_times2)

    # Ecuaciones con términos de acoplamiento
    du[1] = -x1/τ1 + I1 + ε * E2
    du[2] = -x2/τ2 + I2 + ε * E1
end

params = (
    τ1 = 3.0,  # Tiempo de respuesta más rápido para x1
    τ2 = 2.0,  # Tiempo de respuesta mucho más lento para x2
    I1 = 2.0,  # Mayor intensidad positiva para x1
    I2 = -3.0,  # Mayor intensidad negativa para x2
    ε = 0.8,  # Acoplamiento débil
    firing_times1 = [5, 15, 25, 35, 45],  # Más frecuentes para x1
    firing_times2 = [10, 20, 30, 40, 50],  # Menos frecuentes para x2
    std_g = 0.1
)

initial_conditions = [[0.5, 0.6], [1.0, 1.0], [-0.5, -0.6], [0.0, 0.0], [7.0,0.0], [2.0,-11.0],[1.0,4.0]]
tspan = (0.0, 50.0)

# Gráfica de soluciones integrales solo para la primera condición inicial
fig_time = plot(layout=(2, 1), size=(800, 600))
u0 = initial_conditions[1]  # Usar solo la primera condición inicial
prob = ODEProblem(osciladores!, u0, tspan, params)
sol = solve(prob, reltol=1e-6, abstol=1e-6)
plot!(fig_time[1], sol.t, sol[1,:], label="", xlabel=L"t", ylabel=L"x_1(t)", color=:blue)
plot!(fig_time[2], sol.t, sol[2,:], label="", xlabel=L"t", ylabel=L"x_2(t)", color=:blue)

# Campo de fases con múltiples condiciones iniciales
p_phase = plot(title="Retrato Fase", xlabel=L"x_1", ylabel=L"x_2", legend=false)
colors = [:blue, :green, :orange, :purple, :red, :yellow, :pink]  # Colores para diferentes trayectorias


for (i, u0) in enumerate(initial_conditions)
    prob = ODEProblem(osciladores!, u0, tspan, params)
    sol = solve(prob, reltol=1e-6, abstol=1e-6)
    plot!(p_phase, sol[1,:], sol[2,:], label="", color=colors[i])
    scatter!(p_phase, [u0[1]], [u0[2]], label="", color=colors[i], markersize=5)
end

# Mostrar gráficas
display(fig_time)
display(p_phase)