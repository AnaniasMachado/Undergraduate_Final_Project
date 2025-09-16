using CSV
using DataFrames
using MAT
using Base.GC
using Statistics

include("../types.jl")
include("../utility.jl")

problem_folder = "P123_P1Sym"

matrices_folder = "./instances/rectangular_sparse"
m_values = vcat([120 * i for i in 1:4], [600 * i for i in 1:10])

results_folder = "results/problem_$problem_folder"

df = DataFrame()

A_norm_0_list = []
A_norm_1_list = []
AMP_norm_0_list = []
AMP_norm_1_list = []

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
        AMP = pinv(A)

        A_norm_0 = matrix_norm_0(A)
        A_norm_1 = norm(A, 1)
        AMP_norm_0 = matrix_norm_0(AMP)
        AMP_norm_1 = norm(AMP, 1)

        push!(A_norm_0_list, A_norm_0)
        push!(A_norm_1_list, A_norm_1)
        push!(AMP_norm_0_list, AMP_norm_0)
        push!(AMP_norm_1_list, AMP_norm_1)

        GC.gc()

        if idx == max_idx
            A_norm_0_mean = mean(A_norm_0_list)
            A_norm_1_mean = mean(A_norm_1_list)
            AMP_norm_0_mean = mean(AMP_norm_0_list)
            AMP_norm_1_mean = mean(AMP_norm_1_list)

            result = DataFrame(
                m = [m],
                n = [n],
                r = [r],
                d = [d],
                A_norm_0_mean = [A_norm_0_mean],
                A_norm_1_mean = [A_norm_1_mean],
                AMP_norm_0_mean = [AMP_norm_0_mean],
                AMP_norm_1_mean = [AMP_norm_1_mean]
            )

            append!(df, result)

            empty!(A_norm_0_list)
            empty!(A_norm_1_list)
            empty!(AMP_norm_0_list)
            empty!(AMP_norm_1_list)

            GC.gc()
        end
    end
end

results_filename = "data_inst_$(problem_folder).csv"
results_filepath = joinpath(results_folder, results_filename)
CSV.write(results_filepath, df)