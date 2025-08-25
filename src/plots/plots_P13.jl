using CSV
using DataFrames
using Plots
using LaTeXStrings

problem = "P13"

gurobi_filepath = "results/problem_$(problem)/results_$(problem)_Gurobi_PLS.csv"
drs_opt_eps_filepath = "results/problem_$(problem)/results_$(problem)_DRS_Opt_Eps.csv"
drs_opt_r0_filepath = "results/problem_$(problem)/results_$(problem)_DRS_Opt_r0.csv"
drs_fp_eps_filepath = "results/problem_$(problem)/results_$(problem)_DRS_FP_Eps.csv"
drs_fp_r0_filepath = "results/problem_$(problem)/results_$(problem)_DRS_FP_r0.csv"

df_grb = CSV.read(gurobi_filepath, DataFrame)
df_drs_opt_eps = CSV.read(drs_opt_eps_filepath, DataFrame)
df_drs_opt_r0 = CSV.read(drs_opt_r0_filepath, DataFrame)
df_drs_fp_eps = CSV.read(drs_fp_eps_filepath, DataFrame)
df_drs_fp_r0 = CSV.read(drs_fp_r0_filepath, DataFrame)

df_grb = filter(row -> row[:time_mean] != -1, df_grb)
df_drs_opt_eps = filter(row -> row[:time_mean] != -1, df_drs_opt_eps)
df_drs_opt_r0 = filter(row -> row[:time_mean] != -1, df_drs_opt_r0)
df_drs_fp_eps = filter(row -> row[:time_mean] != -1, df_drs_fp_eps)
df_drs_fp_r0 = filter(row -> row[:time_mean] != -1, df_drs_fp_r0)

sort!(df_grb, :m)
sort!(df_drs_opt_eps, :m)
sort!(df_drs_opt_r0, :m)
sort!(df_drs_fp_eps, :m)
sort!(df_drs_fp_r0, :m)

plot(df_grb.m, df_grb.time_mean, label=L"\textrm{Gurobi}", xlabel=L"m", ylabel=L"\textrm{time(s)}")

plot!(df_drs_opt_eps.m, df_drs_opt_eps.time_mean, label=L"\textrm{DRS}_{\textrm{Opt}}^{\epsilon}")

plot!(df_drs_opt_r0.m, df_drs_opt_r0.time_mean, label=L"\textrm{DRS}_{\textrm{Opt}}^{\epsilon(r^0)}")

plot!(df_drs_fp_eps.m, df_drs_fp_eps.time_mean, label=L"\textrm{DRS}_{\textrm{FP}}^{\epsilon}")

plot!(df_drs_fp_r0.m, df_drs_fp_r0.time_mean, label=L"\textrm{DRS}_{\textrm{FP}}^{\epsilon(r^0)}")

# display(current())

savefig("Plots/problem_$(problem)/problem_$(problem)_time_plot.png")

plot(df_grb.m, df_grb.norm_0_ratio_mean, label=L"\textrm{Gurobi}", xlabel=L"m", ylabel=L"||\textrm{H}||_0/||\textrm{AMP}||_0")

plot!(df_drs_opt_eps.m, df_drs_opt_eps.norm_0_ratio_mean, label=L"\textrm{DRS}_{\textrm{Opt}}^{\epsilon}")

plot!(df_drs_opt_r0.m, df_drs_opt_r0.norm_0_ratio_mean, label=L"\textrm{DRS}_{\textrm{Opt}}^{\epsilon(r^0)}")

plot!(df_drs_fp_eps.m, df_drs_fp_eps.norm_0_ratio_mean, label=L"\textrm{DRS}_{\textrm{FP}}^{\epsilon}")

plot!(df_drs_fp_r0.m, df_drs_fp_r0.norm_0_ratio_mean, label=L"\textrm{DRS}_{\textrm{FP}}^{\epsilon(r^0)}")

# display(current())

savefig("Plots/problem_$(problem)/problem_$(problem)_norm_0_ratio_plot.png")

plot(df_grb.m, df_grb.norm_1_ratio_mean, label=L"\textrm{Gurobi}", xlabel=L"m", ylabel=L"||\textrm{H}||_1/||\textrm{AMP}||_1")

plot!(df_drs_opt_eps.m, df_drs_opt_eps.norm_1_ratio_mean, label=L"\textrm{DRS}_{\textrm{Opt}}^{\epsilon}")

plot!(df_drs_opt_r0.m, df_drs_opt_r0.norm_1_ratio_mean, label=L"\textrm{DRS}_{\textrm{Opt}}^{\epsilon(r^0)}")

plot!(df_drs_fp_eps.m, df_drs_fp_eps.norm_1_ratio_mean, label=L"\textrm{DRS}_{\textrm{FP}}^{\epsilon}")

plot!(df_drs_fp_r0.m, df_drs_fp_r0.norm_1_ratio_mean, label=L"\textrm{DRS}_{\textrm{FP}}^{\epsilon(r^0)}")

# display(current())

savefig("Plots/problem_$(problem)/problem_$(problem)_norm_1_ratio_plot.png")