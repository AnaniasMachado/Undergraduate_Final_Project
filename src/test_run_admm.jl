using MAT

include("utility.jl")
include("types.jl")
include("./methods/admm.jl")

m = 100
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

epsilon = 10^(-4)
eps_opt = epsilon
eps_abs = epsilon
eps_rel = 10^(-4)
fixed_tol = true

time_limit = 1200

run_time = @elapsed begin
    H = admm_p123(A, rho, eps_abs, eps_rel, fixed_tol, eps_opt, time_limit)
end
H_norm_0 = matrix_norm_0(H)
H_norm_1 = norm(H, 1)

println("Fixed tol: $fixed_tol")
println("Eps opt: $eps_opt")
println("Eps abs: $eps_abs")
println("Eps rel: $eps_rel")
println("ADMM time: $run_time")
println("ADMM norm 1: $H_norm_1")
println("ADMM norm 0: $H_norm_0")
println("ADMM bound ratio: $(H_norm_0 / (r^2 + (m-r)*n))")