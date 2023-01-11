
module Activity_Binary

function margules(x1, coeff_arr)
    x1_arr = [x1, 1 - x1]

    lnRho = Vector{Float64}(undef, 2)
    for (i, j) in (1:2, 2:-1:1)
        println(i, j)
        lnRho[i] = (coeff_arr[i] + 2 * (coeff_arr[j] - coeff_arr[i]) * x1_arr[i]) * x1_arr[j]^2
    end

    return exp.(lnRho)
end

function van_laar(x1, coeff_arr)

    x1_arr = [x1, 1 - x1]
    A_denom = sum(coeff_arr .* x1_arr)

    lnRho = Vector{Float64}(undef, 2)
    for (i, j) in (1:2, 2:-1:1)
        lnRho[i] = coeff_arr[i] * ((coeff_arr[j] * x1_arr[j]) / (A_denom))^2
    end

    return exp.(lnRho)
end

function wilson(x1, lam_coeff_arr, v_l_arr, T)
    x1_arr = [x1, 1 - x1]

    coeff_arr = Vector{Float64}(undef, 2)
    lnRho = Vector{Float64}(undef, 2)
    
    R=8.314462618
    for (i, j) in (1:2, 2:-1:1)
        coeff_arr[i] = (v_l_arr[j] / v_l_arr[i]) * exp(-lam_coeff_arr[i] / (R * T))
    end

    println(coeff_arr)

    calc_bulk_term = [1, -1] .* (coeff_arr[1] / (x1_arr[1] + coeff_arr[1] * x1_arr[2]) - coeff_arr[2] / (x1_arr[2] + coeff_arr[2] * x1_arr[1]))
    println(calc_bulk_term)

    for (i, j) in (1:2, 2:-1:1)
        lnRho[i] = -log(x1_arr[i] + coeff_arr[i] * x1_arr[j]) + x1_arr[j] * calc_bulk_term[i]
    end

    return exp.(lnRho)

end


end
