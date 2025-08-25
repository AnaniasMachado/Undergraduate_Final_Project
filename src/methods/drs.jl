using LinearAlgebra

function soft_thresholding_matrix(X::Matrix{Float64}, lambda::Float64)
    return sign.(X) .* max.(abs.(X) .- lambda, 0)
end

function get_proj_data(A::Matrix{Float64}, problem::String)
    m, n = size(A)
    if problem == "P1"
        AMP = pinv(A)
        R = A
        S = A
        T = A
        RMP = AMP
        SMP = AMP
        T_factor = RMP * T * SMP
        proj_data = DRSProjDataSimple(R, S, T, RMP, SMP, T_factor, AMP)
        return proj_data
    elseif problem == "P13"
        R = A' * A
        S = I(m)
        T = A'
        RMP = pinv(R)
        SMP = I(m)
        T_factor = RMP * T
        AMP = RMP * A'
        proj_data = DRSProjDataSimple(R, S, T, RMP, SMP, T_factor, AMP)
        return proj_data
    elseif problem == "P14"
        R = I(n)
        S = A * A'
        T = A'
        RMP = I(n)
        SMP = pinv(S)
        T_factor = T * SMP
        AMP = A' * SMP
        proj_data = DRSProjDataSimple(R, S, T, RMP, SMP, T_factor, AMP)
        return proj_data
    elseif (problem in ["P134", "P123", "P124"]) 
        AMP = pinv(A)
        proj_data = DRSProjData(AMP)
        return proj_data
    end
end

function projection(A::Matrix{Float64}, X::Matrix{Float64}, proj_data::Union{DRSProjDataSimple, DRSProjData}, problem::String)
    if problem in ["P123", "P124"]
        Z = proj_data.AMP * A * X * A * proj_data.AMP
        Y = proj_data.AMP - Z + X * A * proj_data.AMP
        return Y
    elseif problem == "P134"
        Z = proj_data.AMP * A * X * A * proj_data.AMP
        Y = X - proj_data.AMP * A * X + proj_data.AMP - X * A * proj_data.AMP + Z
        return Y
    else
        Y = X - proj_data.RMP * proj_data.R * X * proj_data.S * proj_data.SMP + proj_data.T_factor
        return Y
    end
end

function primal_residual(A::Matrix, V::Matrix, proj_data::Union{DRSProjDataSimple, DRSProjData}, problem::String)
    if problem == "P1"
        return norm(A * V * A - A)
    elseif problem == "P13"
        return norm(A' * A * V - A')
    elseif problem == "P14"
        return norm(V * A * A' - A')
    elseif problem in ["P123", "P124"]
        return norm(A' * V' * A' + V * A * proj_data.AMP - A' - V)
    elseif problem == "P134"
        return norm(A' * A * V + V * A * A' - 2 * A')
    end
end

function primal_residual_matrix(A::Matrix, V::Matrix, proj_data::Union{DRSProjDataSimple, DRSProjData}, problem::String)
    if problem == "P1"
        return A * V * A - A
    elseif problem == "P13"
        return A' * A * V - A'
    elseif problem == "P14"
        return V * A * A' - A'
    elseif problem in ["P123", "P124"]
        return A' * V' * A' + V * A * proj_data.AMP - A' - V
    elseif problem == "P134"
        return A' * A * V + V * A * A' - 2 * A'
    end
end

function dual_variable(A::Matrix{Float64}, X::Matrix{Float64}, V::Matrix{Float64}, lambda::Float64, proj_data::Union{DRSProjDataSimple, DRSProjData}, problem::String)
    if problem in ["P123", "P124"]
        B = (1 / lambda) * (X - V)
        L = proj_data.AMP * proj_data.AMP' * B * A * proj_data.AMP
        G = B * A * proj_data.AMP - B
        return L, G
    elseif problem == "P134"
        B = (1 / lambda) * (X - V)
        L = proj_data.AMP * proj_data.AMP' * B
        G = (B - proj_data.AMP * A * B) * proj_data.AMP' * proj_data.AMP
        return L, G
    else
        L = proj_data.RMP' * (1 / lambda) * (X - V) * proj_data.SMP'
        return L
    end
end

function dual_residual(A::Matrix{Float64}, X::Matrix{Float64}, V::Matrix{Float64}, lambda::Float64, proj_data::Union{DRSProjDataSimple, DRSProjData}, problem::String)
    if problem in ["P123", "P124"]
        Lambda, Gamma = dual_variable(A, X, V, lambda, proj_data, problem)
        return norm((1 / lambda) * (V - X) + A' * A * Lambda + Gamma * A * proj_data.AMP - Gamma)
    elseif problem == "P134"
        Lambda, Gamma = dual_variable(A, X, V, lambda, proj_data, problem)
        return norm((1 / lambda) * (V - X) + A' * A * Lambda + Gamma * A * A')
    else
        Lambda = dual_variable(A, X, V, lambda, proj_data, problem)
        return norm((1 / lambda) * (V - X) + proj_data.R' * Lambda * proj_data.S')
    end
end

function dual_residual_matrix(A::Matrix{Float64}, X::Matrix{Float64}, V::Matrix{Float64}, lambda::Float64, proj_data::Union{DRSProjDataSimple, DRSProjData}, problem::String)
    if problem in ["P123", "P124"]
        Lambda, Gamma = dual_variable(A, X, V, lambda, proj_data, problem)
        return (1 / lambda) * (V - X) + A' * A * Lambda + Gamma * A * proj_data.AMP - Gamma
    elseif problem == "P134"
        Lambda, Gamma = dual_variable(A, X, V, lambda, proj_data, problem)
        return (1 / lambda) * (V - X) + A' * A * Lambda + Gamma * A * A'
    else
        Lambda = dual_variable(A, X, V, lambda, proj_data, problem)
        return (1 / lambda) * (V - X) + proj_data.R' * Lambda * proj_data.S'
    end
end

function drs(A::Matrix{Float64}, lambda::Float64, eps_abs::Float64, eps_rel::Float64, problem::String, fixed_tol::Bool, eps_opt::Float64, stop_crit::String, time_limit::Int64)
    start_time = time()
    if problem == "P124"
        A = Matrix(A')
    end
    # Initial data
    m, n = size(A)
    Xh = zeros(n, m)
    X = zeros(n, m)
    # Projection data
    proj_data = get_proj_data(A , problem)
    V = proj_data.AMP
    eps_tol = 10^(-5)
    k = -1
    while true
        k += 1
        Xh = soft_thresholding_matrix(V, lambda)
        Vh = 2 * Xh - V
        X = projection(A, Vh, proj_data, problem)
        if (k == 0) && (fixed_tol == false) && (stop_crit == "Opt")
            initial_pri_res = primal_residual_matrix(A, Xh, proj_data, problem)
            initial_dual_res = dual_res = dual_residual_matrix(A, Xh, V, lambda, proj_data, problem)
            r0 = norm(hcat(initial_pri_res, initial_dual_res))
            eps_tol = eps_abs + eps_rel * r0
        elseif (k == 0) && (fixed_tol == false) && (stop_crit == "Fixed_Point")
            r0 = norm(X - Xh)
            eps_tol = eps_abs + eps_rel * r0
        end
        V += X - Xh
        if fixed_tol && (stop_crit == "Opt")
            pri_res = primal_residual(A, Xh, proj_data, problem)
            dual_res = dual_residual(A, Xh, V, lambda, proj_data, problem)
            if (pri_res <= eps_opt) && (dual_res <= eps_opt)
                break
            end
        elseif !fixed_tol && (stop_crit == "Opt")
            pri_res = primal_residual_matrix(A, Xh, proj_data, problem)
            dual_res = dual_residual_matrix(A, Xh, V, lambda, proj_data, problem)
            res = hcat(pri_res, dual_res)
            if norm(res) <= eps_tol
                break
            end
        elseif fixed_tol && (stop_crit == "Fixed_Point")
            res = norm(X - Xh)
            if res <= eps_opt
                break
            end
        elseif !fixed_tol && (stop_crit == "Fixed_Point")
            res = norm(X - Xh)
            if res <= eps_tol
                break
            end
        end
        # Checks time limit
        elapsed_time = time() - start_time
        if elapsed_time > time_limit
            println("TimeLimit: DRS exceed time limit to solve the problem.")
            return "-", k
        end
    end
    if problem == "P124"
        return Matrix(X'), k
    else
        return X, k
    end
end