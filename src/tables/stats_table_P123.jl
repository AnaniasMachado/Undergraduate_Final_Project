using CSV
using DataFrames
using MAT
using Statistics

include("../utility.jl")

methods = ["Gurobi", "DRS", "DRS_Boyd", "DRS_FP"]

problems = ["P123"]
problem = problems[1]

matrices_folder = "./instances/rectangular_dense_01"

gurobi_solutions_folder = "./solutions/problem_$problem/Gurobi_Cal"

results_folder = "results/problem_$problem"

m_values = [60+i*20 for i in 0:27]
d_values = [100]

df = DataFrame()

for method in methods
    solutions_folder = "./solutions/problem_$problem/$method"

    bound_ratio_list = []
    norm_0_ratio_list = []
    norm_1_ratio_list = []
    rank_ratio_list = []
    m_max = 600
    for m in m_values
        n = Int(m / 2)
        r = Int(m / 4)
        for d in d_values
            for idx in 1:5
                try
                    if method != "Gurobi"
                        sol_name = "problem_$(problem)_m_$(m)_n_$(n)_d_$(d)_idx_5"
                        sol_path = joinpath(solutions_folder, sol_name)
                        sol_data = matread(sol_path)
                    else
                        grb_name = "problem_P123_m_$(m)_n_$(n)_d_$(d)_idx_5"
                        grb_path = joinpath(gurobi_solutions_folder, grb_name)
                        grb_data = matread(grb_path)
                    end
                catch
                    println("Caught error. m = $(m)")
                    m_max = min(m, m_max)
                    break
                end
                mat_name = "experiment_matrix_m$(m)_n$(n)_r$(r)_d$(d)_idx$(idx).mat"
                mat_path = joinpath(matrices_folder, mat_name)
                mat_data = matread(mat_path)
                A = mat_data["matrix"]
                A = Matrix(A)
                AMP = pinv(A)

                H = 0
                if method == "Gurobi"
                    grb_name = "problem_P123_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                    grb_path = joinpath(gurobi_solutions_folder, grb_name)
                    grb_data = matread(grb_path)
                    H = grb_data["H"]
                    H = Matrix(H)
                else
                    sol_name = "problem_$(problem)_m_$(m)_n_$(n)_d_$(d)_idx_$(idx)"
                    sol_path = joinpath(solutions_folder, sol_name)
                    sol_data = matread(sol_path)
                    H = sol_data["H"]
                    H = Matrix(H)
                end

                H_norm_0 = matrix_norm_0(H)
                H_norm_1 = norm(H, 1)
                H_rank = calculate_rank(H)

                H_div_bound_norm_0 = H_norm_0 / (m * r)
                H_div_AMP_norm_0 = H_norm_0 / matrix_norm_0(AMP)
                H_div_AMP_norm_1 = H_norm_1 / norm(AMP, 1)
                H_rank_ratio = H_rank / r

                H_norm_0 = matrix_norm_0(H)
                H_norm_1 = norm(H, 1)

                H_div_AMP_norm_0 = H_norm_0 / matrix_norm_0(AMP)
                H_div_AMP_norm_1 = H_norm_1 / norm(AMP, 1)

                push!(bound_ratio_list, H_div_bound_norm_0)
                push!(norm_0_ratio_list, H_div_AMP_norm_0)
                push!(norm_1_ratio_list, H_div_AMP_norm_1)
                push!(rank_ratio_list, H_rank_ratio)
            end
        end
    end

    bound_ratio_min = minimum(bound_ratio_list)
    bound_ratio_mean = mean(bound_ratio_list)
    bound_ratio_max = maximum(bound_ratio_list)

    norm_0_ratio_min = minimum(norm_0_ratio_list)
    norm_0_ratio_mean = mean(norm_0_ratio_list)
    norm_0_ratio_max = maximum(norm_0_ratio_list)

    norm_1_ratio_min = minimum(norm_1_ratio_list)
    norm_1_ratio_mean = mean(norm_1_ratio_list)
    norm_1_ratio_max = maximum(norm_1_ratio_list)

    rank_ratio_min = minimum(rank_ratio_list)
    rank_ratio_mean = mean(rank_ratio_list)
    rank_ratio_max = maximum(rank_ratio_list)

    result = DataFrame(
        method = [method],
        bound_ratio_min = [bound_ratio_min],
        bound_ratio_mean = [bound_ratio_mean],
        bound_ratio_max = [bound_ratio_max],
        norm_0_ratio_min = [norm_0_ratio_min],
        norm_0_ratio_mean = [norm_0_ratio_mean],
        norm_0_ratio_max = [norm_0_ratio_max],
        norm_1_ratio_min = [norm_1_ratio_min],
        norm_1_ratio_mean = [norm_1_ratio_mean],
        norm_1_ratio_max = [norm_1_ratio_max],
        rank_ratio_min = [rank_ratio_min],
        rank_ratio_mean = [rank_ratio_mean],
        rank_ratio_max = [rank_ratio_max],
        m_max = [m_max]
    )

    append!(df, result)
end

results_filename = "results_$(problem)_stats_table.csv"
results_filepath = joinpath(results_folder, results_filename)
CSV.write(results_filepath, df)