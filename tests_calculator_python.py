# a simple python calculator to find solutions of Volterra equations
# with given parameters. Uses library scipy.integrate.solve_ivp
# to write tests for the CPP program

from scipy.integrate import solve_ivp
import numpy as np

A = 1.3
B = 0.4
C = 0.1
D = 0.5
n_step = 40
delta_t = 0.001
x0 = 231
y0 = 227

sol = solve_ivp(
    lambda t, Y: [(A - B * Y[1]) * Y[0], (C * Y[0] - D) * Y[1]],
    [0, n_step * delta_t],
    [x0, y0],
)
X_n = sol.y[:, -1][0]
Y_n = sol.y[:, -1][1]
H_n = -D * np.log(X_n) + C * X_n + B * Y_n - A * np.log(Y_n)

print(X_n, Y_n, H_n)
