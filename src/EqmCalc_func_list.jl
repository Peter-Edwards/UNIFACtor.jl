
module EqmCalc

include("Activity_modules/Activity_UNIFAC.jl")
include("Activity_modules/Activity_UNIFACmod.jl")
include("Fugacity_modiles/Fugacity_func_list.jl")
include("Options_Input.jl")

using NLopt
const R = 8.314462

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


function func_activity_select(act_method, activity_method_input, T_k, x_arr)

    if act_method == "unifac"
        act_arr = Activity_UNIFAC.UNIFAC(T_k, activity_method_input, x_arr)
    elseif act_method == "unifacmod"
        act_arr = Activity_UNIFACmod.UNIFACmod(T_k, activity_method_input, x_arr)
    end
    return act_arr
end


function func_fugacity_pure_select(fug_method, fug_method_input, T_k, P)

    if fug_method == "unity"
        fug_val = (1, 1)
    elseif fug_method == "vdw"
        fug_val = Fugacity.VdW_EOS(P, T_k, fug_method_input[1], fug_method_input[2])
    elseif fug_method == "rk"
        fug_val = Fugacity.RK_EOS(P, T_k, fug_method_input[1], fug_method_input[2])
    elseif fug_method == "srk"
        fug_val = Fugacity.SRK_EOS(P, T_k, fug_method_input[1], fug_method_input[2], fug_method_input[3])
    elseif fug_method == "pr"
        fug_val = Fugacity.PR_EOS(P, T_k, fug_method_input[1], fug_method_input[2], fug_method_input[3])

    end
    return fug_val
end

function func_fugacity_mix_select(fug_method, fug_method_input, T_k, P, y_arr, kij)
    if fug_method == "unity"
        fug_arr = (1, ones(Float64, length(y_arr)))
    else
        fug_arr = Fugacity.EOS_Mix_full(y_arr, P, T_k, fug_method_input[1], fug_method_input[2], fug_method_input[3], fug_method, kij)
    end
    return fug_arr
end

function func_poynting_add(P, Psat, T, v_l)

    poynting_factor = Vector{Float64}(undef, length(v_l))
    @. poynting_factor = exp(v_l * (P - Psat) / (R * T))

    return poynting_factor
end

function func_bubble_setup(calc_opts, comp_structs)

    #Activity function methods setup
    if calc_opts.activity_mdl == "unifac"
        activity_input = Vector{Matrix{Int64}}(undef, length(comp_structs))
        for i in eachindex(comp_structs)
            activity_input[i] = comp_structs[i].UNIFAC_comp
        end
    elseif calc_opts.activity_mdl == "unifacmod"
        activity_input = Vector{Matrix{Int64}}(undef, length(comp_structs))
        for i in eachindex(comp_structs)
            activity_input[i] = comp_structs[i].UNIFACmod_comp
        end
    end
    #Fugacity function methods setup
    if calc_opts.fugacity_mdl == "unity"
        fugacity_input = Vector{Vector{Float64}}(undef, 0)
        nothing
    elseif calc_opts.fugacity_mdl == "vdw" || calc_opts.fugacity_mdl == "rk"
        fugacity_input = Vector{Vector{Float64}}(undef, 2)
        P_crit_arr = Vector{Float64}(undef, length(comp_structs))
        T_crit_arr = copy(P_crit_arr)
        for i in eachindex(comp_structs)
            P_crit_arr[i] = comp_structs[i].Pc
            T_crit_arr[i] = comp_structs[i].Tc
        end
        fugacity_input[1] = P_crit_arr
        fugacity_input[2] = T_crit_arr

    elseif calc_opts.fugacity_mdl == "srk" || calc_opts.fugacity_mdl == "pr"
        fugacity_input = Vector{Vector{Float64}}(undef, 3)
        P_crit_arr = Vector{Float64}(undef, length(comp_structs))
        T_crit_arr = copy(P_crit_arr)
        ω_crit_arr = copy(P_crit_arr)
        for i in eachindex(comp_structs)
            P_crit_arr[i] = comp_structs[i].Pc
            T_crit_arr[i] = comp_structs[i].Tc
            ω_crit_arr[i] = comp_structs[i].ω
        end
        fugacity_input[1] = P_crit_arr
        fugacity_input[2] = T_crit_arr
        fugacity_input[3] = ω_crit_arr

    end

    #saturation pressure setup
    sat_pressure_input = Vector{Vector{Float64}}(undef, length(comp_structs))
    for i in eachindex(comp_structs)
        sat_pressure_input[i] = comp_structs[i].Antoine_coeff
    end

    if calc_opts.poynting == true
        l_molar_vol = Vector{Float64}(undef, length(comp_structs))
        for i in eachindex(comp_structs)
            l_molar_vol[i] = comp_structs[i].v_l
        end
    else
        l_molar_vol = zeros(Float64, length(comp_structs))
    end

    return activity_input, fugacity_input, sat_pressure_input, l_molar_vol

end


function func_bubble_PXY(optim_vars, comp_structs, calc_opts, T_k, z1_arr, activity_method_input, fugactiy_method_input, sat_pressure_input, l_molar_vol, bub_or_dew,kij)
    #statement to make sure all compositions are possible (faster than using ineqality constraint)
    if sum(optim_vars[2:end]) <= 1
        #extracting temperature and compositions from optim_vars
        P = optim_vars[1]
        z2_arr = Vector{Float64}(undef, length(optim_vars))
        z2_arr[1:length(optim_vars)-1] .= optim_vars[2:end]
        z2_arr[end] = 1 - sum(optim_vars[2:end])

        #bubble or dew determines if x_arr or y_arr is to be found
        if bub_or_dew == "dew"
            x_arr = z2_arr
            y_arr = z1_arr
        else
            x_arr = z1_arr
            y_arr = z2_arr
        end

        #finding the acitivty values for each component
        act_arr = func_activity_select(calc_opts.activity_mdl, activity_method_input, T_k, x_arr)
        #finding the fugacities of each component
        fug_mix_arr = func_fugacity_mix_select(calc_opts.fugacity_mdl, fugactiy_method_input, T_k, P, y_arr, kij)[2]
        #finding the saturation pressure of each component
        P_sat_arr = Antoine(T_k, sat_pressure_input, true)

        fug_pure_arr = Vector{Float64}(undef, length(comp_structs))
        for i in eachindex(comp_structs)
            fug_pure_arr[i] = func_fugacity_pure_select(calc_opts.fugacity_mdl, [fugactiy_method_input[1][i], fugactiy_method_input[2][i], fugactiy_method_input[3][i]], T_k, P_sat_arr[i])[2]
        end

        if calc_opts.poynting == true
            poynting_factor = func_poynting_add(P, P_sat_arr, T_k, l_molar_vol)
            P_sat_arr = P_sat_arr .* poynting_factor
        end

        #defining error for both bub and dew
        ret_val = Vector{Float64}(undef, length(y_arr))
        if bub_or_dew == "dew"
            @. ret_val = x_arr - y_arr * P * fug_mix_arr / (act_arr * P_sat_arr * fug_pure_arr)
        else
            @. ret_val = y_arr - x_arr * act_arr * P_sat_arr * fug_pure_arr / (P * fug_mix_arr)
        end

        return maximum(abs.(ret_val))
    else
        return 1
    end

end

function func_bubble_TXY(optim_vars, comp_structs, calc_opts, P, z1_arr, activity_method_input, fugactiy_method_input, sat_pressure_input, l_molar_vol, bub_or_dew,kij)
    #statement to make sure all compositions are possible (faster than using ineqality constraint)
    if sum(optim_vars[2:end]) <= 1
        #extracting temperature and compositions from optim_vars
        T_k = optim_vars[1]
        z2_arr = Vector{Float64}(undef, length(optim_vars))
        z2_arr[1:length(optim_vars)-1] .= optim_vars[2:end]
        z2_arr[end] = 1 - sum(optim_vars[2:end])

        #bubble or dew determines if x_arr or y_arr is to be found
        if bub_or_dew == "dew"
            x_arr = z2_arr
            y_arr = z1_arr
        else
            x_arr = z1_arr
            y_arr = z2_arr
        end

        #finding the acitivty values for each component
        act_arr = func_activity_select(calc_opts.activity_mdl, activity_method_input, T_k, x_arr)
        #finding the fugacities of each component
        fug_mix_arr = func_fugacity_mix_select(calc_opts.fugacity_mdl, fugactiy_method_input, T_k, P, y_arr, kij)[2]
        #finding the saturation pressure of each component
        P_sat_arr = Antoine(T_k, sat_pressure_input, true)

        fug_pure_arr = Vector{Float64}(undef, length(comp_structs))
        for i in eachindex(comp_structs)
            fug_pure_arr[i] = func_fugacity_pure_select(calc_opts.fugacity_mdl, [fugactiy_method_input[1][i], fugactiy_method_input[2][i], fugactiy_method_input[3][i]], T_k, P_sat_arr[i])[2]
        end

        if calc_opts.poynting == true
            poynting_factor = func_poynting_add(P, P_sat_arr, T_k, l_molar_vol)
            P_sat_arr = P_sat_arr .* poynting_factor
        end

        #defining error for both bub and dew
        ret_val = Vector{Float64}(undef, length(y_arr))
        if bub_or_dew == "dew"
            @. ret_val = x_arr - y_arr * P * fug_mix_arr / (act_arr * P_sat_arr * fug_pure_arr)
        else
            @. ret_val = y_arr - x_arr * act_arr * P_sat_arr * fug_pure_arr / (P * fug_mix_arr)
        end

        return maximum(abs.(ret_val))
    else
        return 1
    end

end


function func_outer_TXY(comp_structs, calc_opts, x_arr, P, bub_or_dew;kij=0)

    activity_input, fugacity_input, sat_pressure_input, l_molar_vol = func_bubble_setup(calc_opts, comp_structs)


    Tb_pure = Antoine(P, sat_pressure_input, false)
    if any(x_arr .>= (1 - 1e-4)) == true
        return [Tb_pure[x_arr.>=(1-1e-4)][1], x_arr[1:(length(x_arr)-1)]]
    else

        function vle_opt_func(x::Vector, grad::Vector)
            if length(grad) > 0
                nothing
            end
            eval = func_bubble_TXY(x, comp_structs, calc_opts, P, x_arr, activity_input, fugacity_input, sat_pressure_input, l_molar_vol, bub_or_dew,kij)
            return eval
        end

        # function func_inequality_composition(x::Vector, grad::Vector)
        #     if length(grad) > 0
        #         nothing
        #     end
        #     return sum(x[2:end]) - (1.0)
        # end

        lb_T = minimum(Tb_pure) * 0.5
        ub_T = maximum(Tb_pure) * 1.5

        opt = Opt(:GN_ORIG_DIRECT_L, length(comp_structs))
        #opt = Opt(:LN_COBYLA, 3)
        #opt=Opt(:GN_AGS, 3)
        #opt=Opt(:GN_ISRES, 3)
        opt.lower_bounds = append!([lb_T], zeros(Float64, length(comp_structs) - 1))
        opt.upper_bounds = append!([ub_T], ones(Float64, length(comp_structs) - 1))
        opt.min_objective = vle_opt_func

        if length(comp_structs) > 2
            #inequality_constraint!(opt::Opt, func_inequality_composition,0)
            opt.maxtime = 5
            stopval!(opt, 1e-3)
        else
            opt.maxtime = 0.2
            stopval!(opt, 1e-4)
        end

        (minf, minx, ret) = optimize(opt, append!([(lb_T + ub_T) / 2], 0.99 / (length(comp_structs)) * (ones(Float64, length(comp_structs) - 1))))
        numevals = opt.numevals # the number of function evaluations
        println("got $minf at $minx after $numevals iterations (returned $ret)")

        return minx
    end
end

function func_outer_PXY(comp_structs, calc_opts, x_arr, T_k, bub_or_dew;kij=0)

    activity_input, fugacity_input, sat_pressure_input, l_molar_vol = func_bubble_setup(calc_opts, comp_structs)


    P_sat_pure = Antoine(T_k, sat_pressure_input, true)
    println(P_sat_pure)
    if any(x_arr .>= (1 - 1e-4)) == true
        return [P_sat_pure[x_arr.>=(1-1e-4)][1], x_arr[1:(length(x_arr)-1)]]
    else

        function vle_opt_func_PXY(x::Vector, grad::Vector)
            if length(grad) > 0
                nothing
            end
            eval = func_bubble_PXY(x, comp_structs, calc_opts, T_k, x_arr, activity_input, fugacity_input, sat_pressure_input, l_molar_vol, bub_or_dew,kij)
            return eval
        end

        # function func_inequality_composition(x::Vector, grad::Vector)
        #     if length(grad) > 0
        #         nothing
        #     end
        #     return sum(x[2:end]) - (1.0)
        # end

        lb_P = minimum(P_sat_pure) * 0.5
        ub_P = maximum(P_sat_pure) * 1.5

        opt = Opt(:GN_ORIG_DIRECT_L, length(comp_structs))
        #opt = Opt(:LN_COBYLA, 3)
        #opt=Opt(:GN_AGS, 3)
        #opt=Opt(:GN_ISRES, 3)
        opt.lower_bounds = append!([lb_P], zeros(Float64, length(comp_structs) - 1))
        opt.upper_bounds = append!([ub_P], ones(Float64, length(comp_structs) - 1))
        opt.min_objective = vle_opt_func_PXY

        if length(comp_structs) > 2
            #inequality_constraint!(opt::Opt, func_inequality_composition,0)
            opt.maxtime = 5
            stopval!(opt, 1e-3)
        else
            opt.maxtime = 0.2
            stopval!(opt, 1e-4)
        end

        (minf, minx, ret) = optimize(opt, append!([(lb_P + ub_P) / 2], 0.99 / (length(comp_structs)) * (ones(Float64, length(comp_structs) - 1))))
        numevals = opt.numevals # the number of function evaluations
        println("got $minf at $minx after $numevals iterations (returned $ret)")

        return minx
    end
end

function gen_graph(comp_structs, calc_opts, P, n_points, mode)

    #this is only meant for 2 components (so it just cuts down comp structs to the first 2 of there are more
    if length(comp_structs) >= 3
        comp_structs = comp_structs[1:2]
    end

    if mode == "fast"
        x_arr = collect(LinRange(0, 1, n_points))
        T_arr = Vector{Float64}(undef, length(x_arr))
        y_arr = Vector{Float64}(undef, length(x_arr))
        for i in 1:n_points
            ret = func_outer_TXY(comp_structs, calc_opts, [x_arr[i], 1 - x_arr[i]], P, "bubble")
            y_arr[i] = ret[2][1]
            T_arr[i] = ret[1]
        end
        return x_arr, y_arr, T_arr
    elseif mode == "full"
        x_arr = collect(LinRange(0, 1, n_points))
        T_bub = Vector{Float64}(undef, length(x_arr))
        T_dew = copy(T_bub)
        for i in 1:n_points
            T_bub[i] = func_outer_TXY(comp_structs, calc_opts, [x_arr[i], 1 - x_arr[i]], P, "bubble")[1]
            T_dew[i] = func_outer_TXY(comp_structs, calc_opts, [x_arr[i], 1 - x_arr[i]], P, "dew")[1]
        end

        return x_arr, T_bub, T_dew
    else
        println("not a mode, choose: fast or full")
    end

end


end
