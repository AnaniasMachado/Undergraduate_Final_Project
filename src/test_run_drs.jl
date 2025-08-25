using MAT

include("utility.jl")
include("types.jl")
include("./methods/drs.jl")

m = 500
n = Int(m / 2)
r = Int(m / 4)
d = 100
idx = 1

matrix_folder = "./instances/rectangular_dense"
mat_file = "A_m$(m)_n$(n)_r$(r)_d$(d)_idx$(idx).mat"
mat_path = joinpath(matrix_folder, mat_file)
mat_data = matread(mat_path)
A = mat_data["A"]
A = Matrix(A)
AMP = pinv(A)
m, n = size(A)
r = rank(A)

println("m = $m, n = $n, r = $r")

data = DataInst(A, m, n, r, AMP=AMP)
constraints = ["P13"]
problem = "P13"
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

DRS_time = @elapsed begin
    DRS_H, DRS_k = drs(A, lambda, eps_abs, eps_rel, problem, fixed_tol, eps_opt, stop_crit, time_limit)
end
DRS_H_norm_0 = matrix_norm_0(DRS_H)
DRS_H_norm_1 = norm(DRS_H, 1)

println("Stop crit: $stop_crit")
println("Fixed tol: $fixed_tol")
println("Eps opt: $eps_opt")
println("Eps abs: $eps_abs")
println("Eps rel: $eps_rel")
println("DRS time: $DRS_time")
println("DRS k: $DRS_k")
println("DRS norm 1: $DRS_H_norm_1")
println("DRS norm 0: $DRS_H_norm_0")
println("DRS bound ratio: $(DRS_H_norm_0 / (m * r))")