using JuMP
import HiGHS

#=  TODO: Validate Dynamic program using LP equivalent.

=#

# SETS
T = 12*24
t = 1:T

# PARAMETERS
E = 240  # max bettery storage (MWh)
r = 10  # max discharge    # TODO Modify to be every 5 minutes (MWh/5min)
s = 10  # max charge       # TODO Modify to be every 5 minutes (MWh/5min)
y0 = 240  # initial battery charge (MWh)

p = [110, 105, 104, 101,  82,  82, 104, 104, 104, 104, 104, 104,
     106, 104, 104, 104, 104, 103, 104, 104, 104, 103, 103, 103,
     103, 103, 103, 103, 100, 100, 102, 102, 101, 102, 101, 101,
     100, 101, 100, 101, 100, 100, 100, 100, 100, 100, 100, 100,
     100, 100, 103, 103, 103, 103, 100, 104, 104, 104, 106, 106,
      83, 100, 103, 103, 108, 109,  83,  84, 101, 105, 117, 129,
      83,  86, 109, 132, 133, 135, 111, 132, 136, 137, 140, 145,
      95, 130, 130, 147, 149, 151, 145, 160, 175, 190, 190, 157,
     200, 151, 149, 149, 139, 130, 149, 149, 149, 131, 113, 113,
     149, 113, 113, 113, 113, 103, 149, 111,  92,  92, 109, 111,
     110,  92,  91,  89,  89,  89,  90,  90,  89,  89,  88,  88,
      90,  90,  90,  89,  89,  87,  77,  80,  77,  81,  84,  74,
      86,  87,  87,  85,  78,  78,  80,  80,  80,  79,  79,  78,
      82,  82,  82,  81,  81,  79,  82,  81,  68,  68,  68,  24,
      67,  59,  59,  58,  58,  58,  45,  45,  45,  45,  35,  34,
      77,  76,  78,  80,  82,  84,  82,  89,  90,  89,  90,  92,
      90,  92,  92,  92,  95, 108,  92,  94,  98, 100, 110, 110,
      92,  92,  94,  97, 108, 110,  95,  98,  99,  99,  98, 104,
     105, 100, 149, 100, 109, 140, 151, 152, 105, 141, 104, 105,
     151, 109, 151, 105,  99,  98, 103,  96,  96,  95,  95,  94,
     113,  94,  93,  91,  91,  90,  95,  95,  95,  92,  92,  91,
      92, 119, 119, 112,  94,  92, 123, 111,  94,  91,  91,  90,
      94,  92,  90,  88,  88,  87, 109,  93,  89,  87,  87,  85,
      93, 109,  93,  85,  84,  84,  86,  92,  85,  84,  84,  83]

# MODEL
model = Model(HiGHS.Optimizer)

# VARIABLES
@variable(model, 0 <= y[t] <= E, Int)
@variable(model, 0 <= u[t] <= r, Int)
@variable(model, 0 <= v[t] <= s, Int)

# EQUATIONS

# objective
@objective(model, Max, sum(p[ti] * (u[ti] - v[ti]) for ti in t))

# constraints
@constraint(model, storage0, y[1] == y0 - u[1] + v[1])
@constraint(model, storage[ti in 2:T], y[ti] .== y[ti-1] - u[ti] + v[ti])

# SOLVE
optimize!(model)
if !is_solved_and_feasible(model)
    error("Solver did not find an optimal solution")
end

# Show solution
objective_value(model)
for i in value(y)
    println("$i")
end