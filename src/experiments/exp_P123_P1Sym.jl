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
problem_folder = "P123_P1Sym"
problems = ["P123", "P1Sym"]
problem = problems[2]
epsilon = 10^(-5)
eps_abs = epsilon
eps_rel = 10^(-3)
fixed_tol = true
eps_opt = epsilon
time_limit = 7200

# Gurobi parameters
constraints_set_P123 = [["P1", "P123", "P3"], ["P13R", "P123"], ["PLS", "P123"], ["PLSr"]]
constraints_P123 = constraints_set_P123[1]

constraints_set_P1Sym = [["P1", "Sym"], ["P1Sym"]]
constraints_P1Sym = constraints_set_P1Sym[1]

constraints_set = [constraints_P123, constraints_P1Sym]
constraints = constraints_set[2]

# DRS parameters
lambda = 10^(-2)

stop_crits = ["Opt", "Fixed_Point"]
stop_crit = stop_crits[2]

matrices_folder = "./instances/rectangular_sparse"
m_values = vcat([120 * i for i in 1:4], [600 * i for i in 1:10])

results_folder = "results/problem_$problem_folder"

solutions_folder = "./solutions/problem_$problem_folder"

df = DataFrame()

norm_0_list = []
norm_1_list = []
rank_list = []
time_list = []

d = 10
idx = 0
max_idx = 5
min_unsolvable_m = Inf

for m in m_values
    for idx in 1:max_idx
        n = Int(m / 3)
        r = Int(m / 4)
        mat_file = "A_m$(m)_n$(n)_r$(r)_d$(d)_idx$(idx).mat"

        mat_path = joinpath(matrices_folder, mat_file)
        mat_data = matread(mat_path)

        println("Solving for matrix: $mat_path")
        
        A = mat_data["A"]
        A = Matrix(A)

        AMP = 0
        data = 0

        if problem == "P123"
            AMP = pinv(A)
            data = DataInst(A, m, n, r, AMP=AMP)
        elseif problem == "P1Sym"
            A = A' * A
            mA, nA = size(A)
            AMP = pinv(A)
            data = DataInst(A, mA, nA, r, AMP=AMP)
        end

        norm_0 = -1.0
        norm_1 = -1.0
        rank = -1.0
        time = -1.0
        if (m < min_unsolvable_m)
            if method == "Gurobi"
                try
                    time = @elapsed begin
                        H = gurobi_solver(data, constraints, eps_opt, time_limit)
                    end
                    norm_0 = matrix_norm_0(H)
                    norm_1 = norm(H, 1)
                    rank = calculate_rank(H)

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
                    norm_0 = matrix_norm_0(H)
                    norm_1 = norm(H, 1)
                    rank = calculate_rank(H)

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
                    norm_0 = matrix_norm_0(H)
                    norm_1 = norm(H, 1)
                    rank = calculate_rank(H)

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

        push!(norm_0_list, norm_0)
        push!(norm_1_list, norm_1)
        push!(rank_list, rank)
        push!(time_list, time)

        GC.gc()

        if idx == max_idx
            norm_0_mean = -1.0
            norm_1_mean = -1.0
            rank_mean = -1.0
            time_mean = -1.0

            if !(-1.0 in norm_0_list)
                norm_0_mean = mean(norm_0_list)
                norm_1_mean = mean(norm_1_list)
                rank_mean = mean(rank_list)
                time_mean = mean(time_list)
            end

            result = DataFrame(
                m = [m],
                n = [n],
                r = [r],
                d = [d],
                norm_0_mean = [norm_0_mean],
                norm_1_mean = [norm_1_mean],
                rank_mean = [rank_mean],
                time_mean = [time_mean]
            )

            append!(df, result)

            empty!(norm_0_list)
            empty!(norm_1_list)
            empty!(rank_list)
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