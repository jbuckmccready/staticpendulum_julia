using ModelingToolkit

@variables t x(t) y(t)
@parameters m g k L d b
D = Differential(t) # define an operator for the differentiation w.r.t. time

gravity_force_x = -m * g * sqrt(1 - (x^2 + y^2) / L^2) / L
gravity_force_y = -m * g * sqrt(1 - (x^2 + y^2) / L^2) / L
drag_force_x = -b * D(x)
drag_force_y = -b * D(y)

eqs = [
    D(D(x)) ~ (gravity_force_x + drag_force_x) / m
    D(D(y)) ~ (gravity_force_y + drag_force_y) / m
]

@named pendulum_sys = ODESystem(eqs, t, [x, y], [m, g, b, k, L, d])

using DifferentialEquations: solve
using Plots: plot

u0 = [
    x => 3,
    y => 3,
    D(x) => 0.0,
    D(y) => 0.0
]

params = [
    d => 0.05, # distance
    m => 1.0, # mass
    g => 9.8, # gravity
    b => 0.1, # drag coefficient
    L => 10.0, # pendulum length
    k => 1.0 # magnetic force coefficient
]

prob = ODEProblem(structural_simplify(pendulum_sys), u0, (0.0, 10.0), params)
sol = solve(prob)
plot(sol, vars=(x, y), aspectratio=1, legend=false)