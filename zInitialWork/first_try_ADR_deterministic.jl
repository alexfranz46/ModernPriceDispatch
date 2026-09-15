using JuMP  #  The equivalent of Python's 'from JuMP import *'
import HiGHS
using DataFrames

#=
Use the JuMP tutorials and the already made GAMS model to create deterministic lookahead model
in Julia.
=#

# SETS
T = 24  # final trading period index
B = 10  # final price tranche index
t = 1:T  # trading periods
b = 1:B  # price tranches

# PARAMETERS

# indexed over b
c = [10, 20, 30, 40, 50, 70, 90, 110, 150, 200]  # $/MWh
bmax = [5, 5, 5, 5, 5, 5, 10, 10, 10, 10]  # MWh

# indexed over t
d = [40, 41, 42, 43, 35, 40, 40, 25, 10, 8, 6, 5, 5, 6, 8, 10, 20, 30, 55, 72, 75, 70, 64, 60]  #MWh

# scalar
Qmax = 70  # MWh
E = 30  # MWh
η = 1
r = 15  # MWh
s = 15  # MWh
rho = 10  # MWh
L = 1000  # $/MWh
q0 = 35  # MWh
y0 = 0  # MWh

# MODEL
model = Model(HiGHS.Optimizer)

# VARIABLES
@variable(model, 0 <= Q[t] <= Qmax, Int)
@variable(model, 0 <= q[t, bi in b] <= bmax[bi], Int)
@variable(model, 0 <= u[t] <= r, Int)
@variable(model, 0 <= v[t] <= s, Int)
@variable(model, 0 <= y[t] <= E, Int)
@variable(model, 0 <= z[ti in t] <= d[ti], Int)

# EQUATIONS

# objective
@objective(model, Min, sum(q[ti,bi] * c[bi] for ti in t, bi in b) + L * sum(z))

# constraints
@constraint(model, setQ, Q[t] .== sum(q[t,bi] for bi in b))
@constraint(model, equilibrium, Q[t] + u[t] - v[t] + z[t] .>= d[t])
@constraint(model, rampUp0, Q[1] - q0 <= rho)
@constraint(model, rampUp, (Q[ti] - Q[ti-1] for ti in 2:T) .<= rho)
@constraint(model, storage0, y[1] == y0 - u[1] + η*v[1])
# @constraint(model, storage, y[ti] .== y[ti-1] - u[ti] + η*v[ti] for ti in 2:T)
@constraint(model, storage[ti in 2:T], y[ti] .== y[ti-1] - u[ti] + η*v[ti])

# SOLVE
optimize!(model)
if !is_solved_and_feasible(model)
    error("Solver did not find an optimal solution")
end


# df = DataFrame(q = value(Q), vminusu=value(v)-value(u))


# termination_status(model)
# primal_status(model)
# dual_status(model)
# objective_value(model)
# value(x)
# value(y)
# shadow_price(c1)
# shadow_price(c2)