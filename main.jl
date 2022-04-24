using ModelingToolkit
using Symbolics: scalarize

@variables t pos[1:2](t)
@parameters m g k L d b
D = Differential(t)

norm_squared = sum(scalarize(pos) .^ 2)
gravity_force = scalarize(pos) .* (-m * g * sqrt(1.0 - norm_squared / L^2) / L)
drag_force = -b .* scalarize(D.(pos))
# drag_force = [-b * D(pos[j]) for j in 1:2]
# gravity_force_x = x * -m * g * sqrt(1 - (x^2 + y^2) / L^2) / L
# gravity_force_y = y * -m * g * sqrt(1 - (x^2 + y^2) / L^2) / L
# drag_force_x = -b * D(x)
# drag_force_y = -b * D(y)

eqs = scalarize(D.(D.(pos)) .~ ((gravity_force .+ drag_force) ./ m))
# eqs = [D(D(pos[j])) ~ (gravity_force[j] + drag_force[j]) / m for j in 1:2]
# eqs = [
#     D(D(x)) ~ (gravity_force_x + drag_force_x) / m
#     D(D(y)) ~ (gravity_force_y + drag_force_y) / m
# ]

@named pendulum_model = ODESystem(eqs, t, [pos...], [m, g, b, k, L, d])
pendulum_sys = structural_simplify(pendulum_model)

using DifferentialEquations: solve
using Plots: plot, @animate, gif

u0 = [
    pos[1] => -8,
    pos[2] => 6,
    D(pos[1]) => 5.0,
    D(pos[2]) => 0.5
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
# plot(sol, vars=(pos[1], pos[2]), aspectratio=1, size=(800, 800), xlims=(-10, 10), ylims=(-10, 10), legend=false)
xyz = map(range(tspan[1], tspan[2], 500)) do t
    x = sol(t)[3]
    y = sol(t)[4]
    z = 10 - sqrt(100 - x^2 - y^2)
    (x, y, z)
end
plot(xyz, aspectratio=1, size=(800, 800), xlabel="x", ylabel="y", zlabel="z", xlims=(-10, 10), ylims=(-10, 10), legend=false, camera=(10, 20))

plt = plot([], [], [], aspectratio=1, size=(800, 800), xlims=(-10, 10), ylims=(-10, 10), zlims=(0, 10), legend=false, camera=(10, 20))
anim = @animate for t = range(tspan[1], tspan[2], 200)
    x = sol(t)[3]
    y = sol(t)[4]
    z = 10 - sqrt(100 - x^2 - y^2)
    push!(plt, x, y, z)
end
gif(anim)