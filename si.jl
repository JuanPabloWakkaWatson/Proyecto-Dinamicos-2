using SymPy

# Define the symbols
@syms x y

# Define the system of equations by setting dx/dt and dy/dt to 0
eq1 = 5/2*x - 1/2*y + 2*x^2 + 1/2*y^2
eq2 = -x + 2*y + 4*x*y

# Solve the system of equations
equilibrium_points = solve([eq1, eq2], [x, y])

# Display the equilibrium points
equilibrium_points