
module EqmCalc

function Antoine(T_or_P, A_params, Psat_mode)


    if Psat_mode == true
        Tb = T_or_P

        P_sat = Array{Float64}(undef, 0)
        for i in eachindex(A_params)
            push!(P_sat, 10^(A_params[i][1] - (A_params[i][2] / (A_params[i][3] + Tb))))
        end
        return P_sat * 100e+3
    else
        Psat = T_or_P * 1 / (100e+3)

        Tb = Array{Float64}(undef, 0)
        for i in eachindex(A_params)
            push!(Tb, A_params[i][2] / (A_params[i][1] - log10(Psat)) - A_params[i][3])
        end
        return Tb
    end
end


using NLopt
include("UNIFAC_func_lstV2.jl")
include("Fugacity_func_list.jl")


function bubb_dew_pt(optim_vars, calc_mode, M_lst, xy_arr, P, Tc, Pc, ω, A_coeff)

    #optim_vars are the required temperature (optim_vars[1]) and unknown mole fractions (optim_vars[2])

    T_k = optim_vars[1]
    if calc_mode == "bubble"
        x_arr = xy_arr
        y_arr = [optim_vars[2], 1 - optim_vars[2]]
    elseif calc_mode == "dew"
        x_arr = [optim_vars[2], 1 - optim_vars[2]]
        y_arr = xy_arr
    else
        x_arr = xy_arr
        y_arr = [optim_vars[2], 1 - optim_vars[2]]
    end

    act_arr = UNIFAC.Activity(T_k, M_lst, x_arr)
    fug_arr = Fugacity.EOS_Mix_full(y_arr, P, T_k, Pc, Tc, ω, "SRK", 0)
    sat_pressure = Antoine(T_k, A_coeff, true)

    pure_fug = Vector{Float64}(undef, length(M_lst))
    for i in eachindex(M_lst)
        pure_fug[i] = Fugacity.SRK_EOS(sat_pressure[i], T_k, Pc[i], Tc[i], ω[i])[2]
    end

    eval = Vector{Float64}(undef, length(y_arr))
    if calc_mode == "bubble"
        @. eval = y_arr - x_arr * act_arr * sat_pressure * pure_fug / (P * fug_arr[2])
    elseif calc_mode == "dew"
        @. eval = x_arr - y_arr * P * fug_arr[2] / (act_arr * sat_pressure * pure_fug)
    end


    return sum((eval) .^ 2) / length(eval)
end



function solve_bubble_dew(calc_mode, M_lst, xy_arr, P, Tc, Pc, ω, A_coeff)

        Tb_pure = Antoine(P, A_coeff, false)
    if any(xy_arr .>= (1 - 1e-4)) == true
        return [Tb_pure[xy_arr .>= (1 - 1e-4)][1], xy_arr[1]]
    else

        function vle_opt_func(x::Vector, grad::Vector)
            if length(grad) > 0
                nothing
            end
            eval = bubb_dew_pt(x, calc_mode, M_lst, xy_arr, P, Tc, Pc, ω, A_coeff)
            return eval

        end
        lb_T=minimum(Tb_pure)*0.8
        ub_T=maximum(Tb_pure)*1.2

        opt = Opt(:GN_ORIG_DIRECT_L, 2)
        #opt = Opt(:LN_COBYLA, 2)
        #opt=Opt(:GN_AGS, 2)
        #opt=Opt(:GN_ISRES, 2)
        opt.lower_bounds = [lb_T, 0]
        opt.upper_bounds = [ub_T, 1]
        opt.stopval = 1e-8
        opt.maxtime = 0.05

        opt.min_objective = vle_opt_func


        (minf, minx, ret) = optimize(opt, [(lb_T+ub_T)/2, 0.5])
        numevals = opt.numevals # the number of function evaluations
        #println("got $minf at $minx after $numevals iterations (returned $ret)")

        return minx
    end

end


function gen_graph(calc_mode, M_lst, P, Tc, Pc, ω, A_coeff, n_points)

    x_arr = collect(LinRange(0, 1, n_points))
    T_arr = Vector{Float64}(undef, length(x_arr))
    Y_arr = Vector{Float64}(undef, length(x_arr))
    for i = 1:n_points
        ret = solve_bubble_dew(
            calc_mode,
            M_lst,
            [x_arr[i], 1 - x_arr[i]],
            P,
            Tc,
            Pc,
            ω,
            A_coeff,
        )
        Y_arr[i] = ret[2]
        T_arr[i] = ret[1]
    end
    return T_arr, Y_arr, x_arr
end


end
