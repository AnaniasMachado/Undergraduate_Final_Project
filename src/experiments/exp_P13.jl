using CSV
using DataFrames
using MAT
using Base.GC
using Statistics

include("../types.jl")
include("../utility.jl")
include("../methods/solvers.jl")
include("../methods/solvers_cal.jl")
include("../methods/drs.jl")

methods = ["Gurobi", "Gurobi_Cal", "DRS"]
method = methods[1]

# Mixed parameters
problems = ["P13"]
problem = problems[1]
epsilon = 10^(-5)
eps_abs = epsilon
eps_rel = 10^(-4)
fixed_tol = false
eps_opt = epsilon
time_limit = 7200

# Gurobi parameters
constraints_set = [["P1", "P3"], ["P13R"], ["PLS"]]
constraints = constraints_set[1]

# DRS parameters
lambda = 10^(-2)

stop_crits = ["Opt", "Fixed_Point"]
stop_crit = stop_crits[1]

matrices_folder = "./instances/rectangular_dense"
m_values = [500 * i for i in 1:10]

results_folder = "results/problem_$problem"

solutions_folder = "./solutions/problem_$problem"

df = DataFrame()

bound_ratio_list = []
norm_0_ratio_list = []
norm_1_ratio_list = []
rank_ratio_list = []
time_list = []

d = 100
max_idx = 5
min_unsolvable_m = Inf

for m in m_values
    for idx in 1:max_idx
        n = Int(m / 2)
        r = Int(m / 4)
        mat_file = "A_m$(m)_n$(n)_r$(r)_d$(d)_idx$(idx).mat"

        mat_path = joinpath(matrices_folder, mat_file)
        mat_data = matread(mat_path)

        println("Solving for matrix: $mat_path")
        
        A = mat_data["A"]
        A = Matrix(A)
        AMP = pinv(A)

        data = DataInst(A, m, n, r, AMP=AMP)

        bound_ratio = -1.0
        norm_0_ratio = -1.0
        norm_1_ratio = -1.0
        rank_ratio = -1.0
        time = -1.0
        if (m < min_unsolvable_m)
            if method == "Gurobi"
                try
                    time = @elapsed begin
                        H = gurobi_solver(data, constraints, eps_opt, time_limit)
                    end
                    H_norm_0 = matrix_norm_0(H)
                    H_norm_1 = norm(H, 1)
                    H_rank = calculate_rank(H)

                    bound_ratio = H_norm_0 / (m * r)
                    norm_0_ratio = H_norm_0 / matrix_norm_0(AMP)
                    norm_1_ratio = H_norm_1 / norm(AMP, 1)
                    rank_ratio = H_rank / r

                    problem_label = join(constraints, "_")

                    solution_filename = "Gurobi/problem_$(problem_label)_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                    solution_filepath = joinpath(solutions_folder, solution_filename)
                    matwrite(solution_filepath, Dict("H" => H, "time" => time))
                catch e
                    if isa(e, ErrorException)
                        global min_unsolvable_m = min(m, min_unsolvable_m)
                    else
                        throw(ErrorException("Gurobi failed to solve problem something unexpected.", e))
                    end
                end
            elseif method == "Gurobi_Cal"
                try
                    time = @elapsed begin
                        H = gurobi_solver_cal(data, problem, eps_opt, time_limit)
                    end
                    H_norm_0 = matrix_norm_0(H)
                    H_norm_1 = norm(H, 1)
                    H_rank = calculate_rank(H)

                    bound_ratio = H_norm_0 / (m * r)
                    norm_0_ratio = H_norm_0 / matrix_norm_0(AMP)
                    norm_1_ratio = H_norm_1 / norm(AMP, 1)
                    rank_ratio = H_rank / r

                    solution_filename = "Gurobi_Cal/problem_$(problem)_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                    solution_filepath = joinpath(solutions_folder, solution_filename)
                    matwrite(solution_filepath, Dict("H" => H, "time" => time))
                catch e
                    if isa(e, ErrorException)
                        global min_unsolvable_m = min(m, min_unsolvable_m)
                    else
                        throw(ErrorException("Gurobi failed to solve problem something unexpected.", e))
                    end
                end
            elseif method == "DRS"
                time = @elapsed begin
                    H, k = drs(A, lambda, eps_abs, eps_rel, problem, fixed_tol, eps_opt, stop_crit, time_limit)
                end
                if H == "-"
                    global min_unsolvable_m = min(m, min_unsolvable_m)
                else
                    H_norm_0 = matrix_norm_0(H)
                    H_norm_1 = norm(H, 1)
                    H_rank = calculate_rank(H)

                    bound_ratio = H_norm_0 / (m * r)
                    norm_0_ratio = H_norm_0 / matrix_norm_0(AMP)
                    norm_1_ratio = H_norm_1 / norm(AMP, 1)
                    rank_ratio = H_rank / r

                    if fixed_tol && stop_crit == "Opt"
                        solution_filename = "DRS_Opt_Eps/problem_$(problem)_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                        solution_filepath = joinpath(solutions_folder, solution_filename)
                        matwrite(solution_filepath, Dict("H" => H, "time" => time, "k" => k))
                    elseif !fixed_tol && stop_crit == "Opt"
                        solution_filename = "DRS_Opt_r0/problem_$(problem)_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                        solution_filepath = joinpath(solutions_folder, solution_filename)
                        matwrite(solution_filepath, Dict("H" => H, "time" => time, "k" => k))
                    elseif fixed_tol && stop_crit == "Fixed_Point"
                        solution_filename = "DRS_FP_Eps/problem_$(problem)_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                        solution_filepath = joinpath(solutions_folder, solution_filename)
                        matwrite(solution_filepath, Dict("H" => H, "time" => time, "k" => k))
                    elseif !fixed_tol && stop_crit == "Fixed_Point"
                        solution_filename = "DRS_FP_r0/problem_$(problem)_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                        solution_filepath = joinpath(solutions_folder, solution_filename)
                        matwrite(solution_filepath, Dict("H" => H, "time" => time, "k" => k))
                    end
                end
            else
                throw(ErrorException("Invalid method chose."))
            end
        end

        push!(bound_ratio_list, bound_ratio)
        push!(norm_0_ratio_list, norm_0_ratio)
        push!(norm_1_ratio_list, norm_1_ratio)
        push!(rank_ratio_list, rank_ratio)
        push!(time_list, time)

        GC.gc()

        if idx == max_idx
            bound_ratio_mean = -1.0
            norm_0_ratio_mean = -1.0
            norm_1_ratio_mean = -1.0
            rank_ratio_mean = -1.0
            time_mean = -1.0

            if !(-1.0 in bound_ratio_list)
                bound_ratio_mean = mean(bound_ratio_list)
                norm_0_ratio_mean = mean(norm_0_ratio_list)
                norm_1_ratio_mean = mean(norm_1_ratio_list)
                rank_ratio_mean = mean(rank_ratio_list)
                time_mean = mean(time_list)
            end

            result = DataFrame(
                m = [m],
                n = [n],
                r = [r],
                d = [d],
                bound_ratio_mean = [bound_ratio_mean],
                norm_0_ratio_mean = [norm_0_ratio_mean],
                norm_1_ratio_mean = [norm_1_ratio_mean],
                rank_ratio_mean = [rank_ratio_mean],
                time_mean = [time_mean]
            )

            append!(df, result)

            empty!(bound_ratio_list)
            empty!(norm_0_ratio_list)
            empty!(norm_1_ratio_list)
            empty!(rank_ratio_list)
            empty!(time_list)

            GC.gc()
        end
    end
end

if method == "Gurobi"
    problem_label = join(constraints, "_")
    results_filename = "results_$(problem_label)_Gurobi.csv"
    results_filepath = joinpath(results_folder, results_filename)
    CSV.write(results_filepath, df)
elseif method == "Gurobi_Cal"
    results_filename = "results_$(problem)_Gurobi_Cal.csv"
    results_filepath = joinpath(results_folder, results_filename)
    CSV.write(results_filepath, df)
elseif method == "DRS"
    if fixed_tol && stop_crit == "Opt"
        results_filename = "results_$(problem)_DRS_Opt_Eps.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    elseif !fixed_tol && stop_crit == "Opt"
        results_filename = "results_$(problem)_DRS_Opt_r0.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    elseif fixed_tol && stop_crit == "Fixed_Point"
        results_filename = "results_$(problem)_DRS_FP_Eps.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    else
        results_filename = "results_$(problem)_DRS_FP_r0.csv"
        results_filepath = joinpath(results_folder, results_filename)
        CSV.write(results_filepath, df)
    end
else
    throw(ErrorException("Invalid method chose."))
end