#Pregunta 2

using DifferentialEquations
using Plots
using LaTeXStrings
using NLsolve
using LinearAlgebra

# Define las ecuaciones del sistema
function sistema_ecuaciones!(du, u, p, t)
    x, y = u
    du[1] = 5/2 * x - 1/2 * y + 2 * x^2 + 1/2 * y^2
    du[2] = -x + 2 * y + 4 * x * y
end

# Puntos de equilibrio: resolver du/dt = 0 para x e y
function encontrar_puntos_equilibrio()
    puntos_equilibrio = []


    function f!(F, u)
        x, y = u
        F[1] = 5/2 * x - 1/2 * y + 2 * x^2 + 1/2 * y^2
        F[2] = -x + 2 * y + 4 * x * y
    end

    # Estimaciones iniciales
    initial_guesses = [[0.0, 0.0], [1.0, 1.0], [-1.0, -1.0], [0.0, 3.0], [3.0, 0.0], [-3.0,-3.0]]
    for guess in initial_guesses
        result = nlsolve(f!, guess)
        if converged(result)
            punto_equilibrio = result.zero
            # Asegurarse de que no se añada un punto ya encontrado debido a la precisión numérica
            if all(!isapprox(punto, punto_equilibrio, atol=1e-5) for punto in puntos_equilibrio)
                push!(puntos_equilibrio, punto_equilibrio)
            end
        end
    end

    return puntos_equilibrio
end

# Calcular y mostrar los puntos de equilibrio
puntos_equilibrio = encontrar_puntos_equilibrio()
println("Puntos de equilibrio encontrados: ", puntos_equilibrio)

# Análisis de estabilidad a través de la matriz Jacobiana
function analisis_estabilidad(punto)
    x, y = punto
    J = [5/2 + 4*x  -1/2;
         4*y        2 + 4*x]
    autovalores = eigvals(J)
    return autovalores
end

# Clasificación de la estabilidad de los puntos de equilibrio
for punto in puntos_equilibrio
    autovalores = analisis_estabilidad(punto)
    println("Punto de equilibrio: ", punto, " Autovalores: ", autovalores)
    
    if all(λ -> real(λ) < 0, autovalores)
        println("El punto ", punto, " es estable (nodo atractor o espiral).")
    elseif all(λ -> real(λ) > 0, autovalores)
        println("El punto ", punto, " es inestable (nodo repulsor o espiral).")
    elseif any(λ -> real(λ) == 0, autovalores)
        if all(imag(λ) ≠ 0 for λ in autovalores)
            println("El punto ", punto, " es un centro.")
        else
            println("El punto ", punto, " puede ser un punto silla o requiere análisis no lineal.")
        end
    else
        println("El punto ", punto, " es un punto silla.")
    end
end




#-------------------------------------------------------------------------------------------------------------

#Pregunta 14

# Define los parámetros del sistema
m = 1.0   # masa
b1 = 0.1  # parámetro de resistencia lineal
b2 = 0.1  # parámetro de resistencia no lineal de primer orden
b3 = 0.1  # parámetro de resistencia no lineal de segundo orden
b4 = 0.2  # parámetro de rigidez no lineal de primer orden; diferente de b3*m
b5 = 0.1  # parámetro de masa efectiva no lineal
b6 = 0.1  # parámetro de rigidez no lineal de segundo orden

# Define las ecuaciones del sistema
function sistema_ecuaciones14(du, u, p, t)
    x, v = u
    du[1] = v
    du[2] = -(b1*x + b2*x*v + b3*x^2*v + b4*x^3 + b6*x^3) / (m + b5*x^2)
end

# Condiciones iniciales
x0 = 2 / (-b3/b2 + b4/(m*b2))
v0 = 0.0
u0 = [x0, v0]

# Define el rango de tiempo para la simulación
tspan = (0.0, 50.0)

# Define el problema de ODE con las condiciones iniciales y el rango de tiempo
prob = ODEProblem(sistema_ecuaciones14, u0, tspan)

# Resuelve el problema
sol = solve(prob)

# Grafica el retrato de fase
plot(sol, vars=(1,2), xlabel=L"x", ylabel=L"v", label="Órbita", title="Retrato de Fase")
# Identificar y graficar la condición inicial
scatter!([x0], [v0], color=:red, markersize=5, label="Condición inicial")

# Muestra el gráfico
display(plot)

#------------------------------------------------------------------------------------------------------
#Pregunta 15
#Pregunta 15

# Define las ecuaciones del sistema no lineal
function sistema_ecuaciones15(du, u, p, t)
    x, y = u
    a1, n1, b1, a2, n2, b2 = p
    du[1] = a1 / (1 + y^n1) + b1 * y
    du[2] = a2 / (1 + x^n2) + b2 * x
end

# Parámetros para el primer caso
p1 = (4.0, 3, -1.0, 4.0, 3, -1.0)
p2 = (4.0, 1, -1.0, 4.0, 3, 1.0)

# Encuentra los puntos de equilibrio
function encontrar_puntos_equilibrio(p)
    function f!(F, u)
        x, y = u
        F[1] = p[1] / (1 + y^p[2]) + p[3] * y
        F[2] = p[4] / (1 + x^p[5]) + p[6] * x
    end
    sol = nlsolve(f!, [0.0, 0.0])
    return sol.zero
end

punto_equilibrio = encontrar_puntos_equilibrio(p2)

# Define el rango de tiempo para la simulación
tspan = (0.0, 10.0)

# Condiciones iniciales aleatorias para la visualización del retrato de fase
u0 = [0.5, 0.5]

# Define el problema de ODE con las condiciones iniciales y el rango de tiempo
prob = ODEProblem(sistema_ecuaciones15, u0, tspan, p2)

# Resuelve el problema
sol = solve(prob)

# Grafica el retrato de fase
plot(sol, vars=(1,2), xlabel=L"x", ylabel=L"y", label="Trayectoria", title="Retrato de Fase")
scatter!([punto_equilibrio[1]], [punto_equilibrio[2]], color=:red, markersize=5, label="Punto de equilibrio")

# Muestra el gráfico
display(plot)
#------------------------------------------------------------------------------------------------------------------------------
#Pregunta 16
# Parámetros del sistema
ω0 = 1.0
ϵ_values = [0.45, 0.25, 0.01]

# Define las ecuaciones del sistema
function sistema_ecuaciones16(du, u, p, t)
    x, y = u
    ϵ, ω0 = p
    du[1] = y
    du[2] = -ϵ * abs(y) * y - ω0^2 * x
end

# Condiciones iniciales y rango de tiempo
u0 = [2.0, 0.0]
tspan = (0.0, 50.0)

# Resuelve y grafica para cada valor de ϵ
for ϵ in ϵ_values
    # Define el problema de ODE con las condiciones iniciales, rango de tiempo y parámetros
    prob = ODEProblem(sistema_ecuaciones16, u0, tspan, (ϵ, ω0))

    # Resuelve el problema
    sol = solve(prob)

    # Gráficas
    p1 = plot(sol, vars=(1,2), label="y vs. x", title="Espacio de fase", xlabel="x", ylabel="y")
    p2 = plot(sol.t, getindex.(sol.u, 1), label="x(t)", title="Posición vs. Tiempo", xlabel="t", ylabel="x")
    p3 = plot(sol.t, getindex.(sol.u, 2), label="y(t)", title="Velocidad vs. Tiempo", xlabel="t", ylabel="y")
    p4 = plot(sol.t, [u[1]^2 + u[2]^2 for u in sol.u], label="Energía", title="Energía vs. Tiempo", xlabel="t", ylabel="Energía")

    # Combinar las gráficas en una figura con subgráficas
    combined_plot = plot(p1, p2, p3, p4, layout=(2,2), legend=true)

    # Mostrar las gráficas
    display(combined_plot)
end