using MAT

include("utility.jl")
include("types.jl")
include("./methods/solvers.jl")
include("./methods/solvers_cal.jl")
include("./methods/drs.jl")

m = 40
n = m
r = Int(m / 4)
d = 100
idx = 1

A = gen_random_rank_r_matrix(m, n, r)

# matrix_folder = "./instances/square_dense"
# mat_file = "A_m$(m)_n$(n)_r$(r)_d$(d)_idx$(idx).mat"
# mat_path = joinpath(matrix_folder, mat_file)
# mat_data = matread(mat_path)
# A = mat_data["A"]
# A = Matrix(A)

A = A' * A
AMP = pinv(A)
m, n = size(A)
r = rank(A)

println("m = $m, n = $n, r = $r")

problem = "P1Sym"
data = DataInst(A, m, n, r, AMP=AMP)
constraints = ["P1", "Sym"]
rho = 3.0
lambda = 10^(-2)

epsilon = 10^(-5)
eps_opt = epsilon
eps_abs = epsilon
eps_rel = 10^(-4)
fixed_tol = false

time_limit = 1200
stop_crits = ["Opt", "Fixed_Point"]
stop_crit = stop_crits[1]

GRB_time = @elapsed begin
    GRB_H = gurobi_solver(data, constraints, eps_opt, time_limit)
end
GRB_H_norm_0 = matrix_norm_0(GRB_H)
GRB_H_norm_1 = norm(GRB_H, 1)

# GRB_Cal_time = @elapsed begin
#     GRB_Cal_H = gurobi_solver_cal(data, problem, eps_opt, time_limit)
# end
# GRB_Cal_H_norm_0 = matrix_norm_0(GRB_Cal_H)
# GRB_Cal_H_norm_1 = norm(GRB_Cal_H, 1)

DRS_time = @elapsed begin
    DRS_H, DRS_k = drs(A, lambda, eps_abs, eps_rel, problem, fixed_tol, eps_opt, stop_crit, time_limit)
end
DRS_H_norm_0 = matrix_norm_0(DRS_H)
DRS_H_norm_1 = norm(DRS_H, 1)

println("GRB time: $GRB_time")
println("GRB norm 1: $GRB_H_norm_1")
println("GRB norm 0: $GRB_H_norm_0")

# println("GRB Cal time: $GRB_Cal_time")
# println("GRB Cal norm 1: $GRB_Cal_H_norm_1")
# println("GRB Cal norm 0: $GRB_Cal_H_norm_0")

println("Stop crit: $stop_crit")
println("Fixed tol: $fixed_tol")
println("Eps opt: $eps_opt")
println("Eps abs: $eps_abs")
println("Eps rel: $eps_rel")
println("DRS time: $DRS_time")
println("DRS k: $DRS_k")
println("DRS norm 1: $DRS_H_norm_1")
println("DRS norm 0: $DRS_H_norm_0")