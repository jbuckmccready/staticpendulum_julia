using ModelingToolkit

@variables t x(t) y(t)
@parameters m g k L d b
D = Differential(t)

gravity_force_x = x * -m * g * sqrt(1 - (x^2 + y^2) / L^2) / L
gravity_force_y = y * -m * g * sqrt(1 - (x^2 + y^2) / L^2) / L
drag_force_x = -b * D(x)
drag_force_y = -b * D(y)

eqs = [
    D(D(x)) ~ (gravity_force_x + drag_force_x) / m
    D(D(y)) ~ (gravity_force_y + drag_force_y) / m
]

@named pendulum_model = ODESystem(eqs, t, [x, y], [m, g, b, k, L, d])
pendulum_sys = structural_simplify(pendulum_model)

equations(pendulum_sys)

using DifferentialEquations: solve
using Plots: plot, @animate, gif

u0 = [
    x => 3,
    y => 7,
    D(x) => 5.0,
    D(y) => 0.5
]

params = [
    d => 0.05, # distance
    m => 1.0, # mass
    g => 9.8, # gravity
    b => 0.1, # drag coefficient
    L => 10.0, # pendulum length
    k => 1.0 # magnetic force coefficient
]

tspan = (0.0, 75.0)

prob = ODEProblem(pendulum_sys, u0, tspan, params)
sol = solve(prob)
plot(sol, vars=(x, y), aspectratio=1, size=(800, 800), xlims=(-10, 10), ylims=(-10, 10), legend=false)

# sol(0)
# sol(0)[2]

plt = plot([], [], aspectratio=1, size=(800, 800), xlims=(-10, 10), ylims=(-10, 10), legend=false)
anim = @animate for t = range(tspan[1], tspan[2], 300)
    push!(plt, sol(t)[3], sol(t)[4])
end
gif(anim)