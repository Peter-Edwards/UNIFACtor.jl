

module Activity_UNIFACmod

export UNIFACmod

#the overall activity function
function UNIFACmod(T_k::Float64, M_lst::Vector{Matrix{Int64}}, x_arr::Vector{Float64})

    #changes zero values so the the math is doable. useful for generating PXY graphs and such.
    for i in eachindex(x_arr)
        if x_arr[i] <= 1e-8
            x_arr[i] = 1e-8
        end
    end

    #--residual contribution
    #forming a matrix of group information for mixture and pure substance interraction
    grp_tab::InfoStruct = func_infotab_grpV3(M_lst, x_arr)
    ind_tab::Vector{InfoStruct} = func_infotab_indV3(M_lst)

    #updating each struct with the LnRho line (LnRHO being the natural log of the activity coefficients of a group 
    #in compound i.
    grp_tab = func_lnRho_V2(grp_tab, T_k)
    for i in eachindex(ind_tab)
        ind_tab[i] = func_lnRho_V2(ind_tab[i], T_k)
    end

    #uses the ln() group activity coefficients for the individual components and the group as a whole 
    # to find the residual ln() residual contribution for a molecule.
    idx = Array{Float64}(undef, length(ind_tab))
    for i in eachindex(ind_tab)
        idx[i] = func_lngamV2(ind_tab[i], grp_tab)
    end

    #--combining combinatorial and residual contributions
    gam_arr = exp.(func_comb_mk2(x_arr, M_lst) + idx)

    return gam_arr

end

#finding Ri and Qi for every molecule
function func_qiri_mk3(M_local::Vector{Matrix{Int64}})
    Ri = zeros(Float64, length(M_local))
    Qi = zeros(Float64, length(M_local))

    for i in eachindex(M_local)
        for j in axes(M_local[i], 1)
            Ri[i] += Rk_arr_glob[M_local[i][j, 1]] * M_local[i][j, 2]
            Qi[i] += Qk_arr_glob[M_local[i][j, 1]] * M_local[i][j, 2]
        end
    end

    return (Ri, Qi)
end

#finding theta and phi
function func_thetaphi(Ri::Array, Qi::Array, x_arr::Array)
    phi::Vector{Float64} = Ri ./ sum(Ri .* x_arr)
    theta::Vector{Float64} = Qi ./ sum(Qi .* x_arr)

    phi_mod::Vector{Float64} = Ri .^ (3 / 4) ./ sum(Ri .^ (3 / 4) .* x_arr)

    return (phi, theta, phi_mod)
end

#finding the overall combinatorial contributions
function func_comb_mk2(x_arr, M_local)
    #(Ri, Qi) = func_qiri_mk2(M_local,Qk_arr_glob,Rk_arr_glob)
    (Ri, Qi) = func_qiri_mk3(M_local)
    (phi, theta, phi_mod) = func_thetaphi(Ri, Qi, x_arr)

    ln_gamc = Vector{Float64}(undef, length(Ri))

    @. ln_gamc = 1 - phi_mod + log(phi_mod) - 5 * Qi * (1 - phi / theta + log(phi / theta))
    return ln_gamc
end

#the func_idxarr function is just used to generate the group indexes for f_idxarr
# function func_idxarr()
#     size_arr=[4;5;2;3;3;1;1;1;2;1;2;1;3;4;3;2;1;3;2;1;3;3;1;1;1;3;1;1;2;1;1;1;1;2;1;1;1;1;2;3;1;3;3;1;1;4;3;3;2;3;5;3;2;3;2;1;1;1;2;1;1;2;1]
#     idx_arr=Array{Float64 }(undef,0)
#     ind=1
#     for i in 1:length(size_arr)
#         append!(idx_arr,ones(spush!(ind.index,1)arr[i])*ind)
#         ind=ind+1
#     end
#     return idx_arr
# end


#find area fraction for each pure molecule and for the mixture
function func_arrfrac(Qk, M_lst, x_arr)
    #M_lst=M_lst.*x_arr'
    M_lst = M_lst .* x_arr
    M_ex = sum(M_lst, dims=1)
    M_ex = M_ex / sum(M_ex)
    Ar_frac = M_ex .* Qk'
    Ar_denom = sum(Ar_frac)
    Ar_frac = Ar_frac ./ Ar_denom

    return Ar_frac


end


#the struct for infotab
mutable struct InfoStruct

    index::Array{Int64}
    Quan::Array{Int64}
    Group::Array{Int64}
    AreaFrac::Array{Float64}
    Qk::Array{Float64}
    LnRho::Array{Float64}

end

function func_arrfrac_indV2(M_local_comp)
    frac_comp = Vector{Float64}(undef, size(M_local_comp, 1))
    frac_comp =
        M_local_comp[:, 2] ./ sum(M_local_comp[:, 2]) .* Qk_arr_glob[M_local_comp[:, 1]]

    return frac_comp / sum(frac_comp)
end


function func_infotab_indV3(M_local)
    tab_out = Vector{InfoStruct}(undef, length(M_local))

    for i in eachindex(M_local)
        struc_frac = M_local[i][:, 1]
        Ar_frac = func_arrfrac_indV2(M_local[i])

        tab_out[i] = InfoStruct(
            struc_frac,
            M_local[i][:, 2],
            Idxarr_glob[struc_frac],
            Ar_frac,
            Qk_arr_glob[struc_frac],
            [],
        )

    end
    return tab_out
end


function func_infotab_grpV3(M_local, x_local)

    join_M_local = reduce(vcat, M_local)
    struct_frac = unique(join_M_local[:, 1])
    Quan_struc = zeros(Float64, length(struct_frac))
    for i in axes(join_M_local, 1)
        for j in axes(struct_frac)
            if join_M_local[:, 1][i] == struct_frac[j]
                Quan_struc[j] += join_M_local[j, 2]
            end
        end
    end


    prealloc_arrfrac = zeros(Float64, length(struct_frac))
    for i in eachindex(M_local)
        for j in axes(M_local[i], 1)
            prealloc_arrfrac[M_local[i][j, 1].==struct_frac] .+=
                M_local[i][j, 2] * x_local[i]
        end
    end
    Xm_arr::Vector{Float64} = prealloc_arrfrac ./ sum(prealloc_arrfrac)
    Ar_frac::Vector{Float64} =
        Qk_arr_glob[struct_frac] .* Xm_arr ./ sum(Qk_arr_glob[struct_frac] .* Xm_arr)

    return InfoStruct(
        struct_frac,
        Quan_struc,
        Idxarr_glob[struct_frac],
        Ar_frac,
        Qk_arr_glob[struct_frac],
        [],
    )
end


function func_lnRho_V2(info_in::InfoStruct, T::Float64)

    num_comp::Int64 = length(info_in.index)
    psi_mat = Matrix{Float64}(undef, num_comp, num_comp)
    for j = 1:num_comp
        for i = 1:num_comp
            if j == i
                psi_mat[i, j] = one(Float64)
            else
                A_nm_var = Anm[info_in.Group[i]+1, info_in.Group[j]+1]
                B_nm_var = Bnm[info_in.Group[i]+1, info_in.Group[j]+1]
                C_nm_var = Cnm[info_in.Group[i]+1, info_in.Group[j]+1]
                psi_mat[i, j] = exp(-(A_nm_var / T + B_nm_var + C_nm_var * T))
            end

        end
    end


    denom_sums = Vector{Float64}(undef, length(info_in.index))
    for i in eachindex(denom_sums)
        denom_sums[i] = sum(info_in.AreaFrac .* view(psi_mat, :, i))
    end

    LnRHO = Vector{Float64}(undef, length(info_in.index))
    pre_multi = Vector{Float64}(undef, length(info_in.AreaFrac))
    pre_multi = info_in.AreaFrac ./ denom_sums

    for i in eachindex(info_in.index)
        LnRHO[i] =
            info_in.Qk[i] *
            (1 - log(denom_sums[i]) - sum(pre_multi .* view(psi_mat, i, :)))
    end
    info_in.LnRho = LnRHO
    return info_in
end


#combines lnRhos for each group from the pure component and the mixture. does this one molecule at a time
function func_lngamV2(info_ind::InfoStruct, info_grp::InfoStruct)

    n = info_ind.index .∈ [info_grp.index]
    idx_grp = indexin(info_ind.index[n], info_grp.index)
    idx_ind = indexin(info_ind.index[n], info_ind.index)

    lngamR = zero(Float64)
    for i in eachindex(idx_ind)
        lngamR +=
            (info_grp.LnRho[idx_grp[i]] - info_ind.LnRho[idx_ind[i]]) *
            info_ind.Quan[idx_ind[i]]
    end

    return lngamR
end



#parameter data

Anm = [
    0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 16.0 17.0 18.0 19.0 20.0 21.0 22.0 23.0 24.0 25.0 26.0 27.0 28.0 29.0 30.0 31.0 32.0 33.0 34.0 35.0 36.0 37.0 38.0 39.0 40.0 41.0 42.0 43.0 44.0 45.0 46.0 47.0 48.0 49.0 52.0 53.0 55.0 56.0 61.0 84.0 85.0 87.0 89.0 90.0 91.0 93.0 98.0 99.0
    1.0 0.0 189.66 114.2 7.339 2777.0 2409.4 1391.3 1381.0 433.6 875.85 98.65601 508.4 233.1 -164.04 350.58 -175.7 958.74 1802.3 593.07 1182.2 401.0 -233.66 -653.74 267.51 -1385.0 2345.0 2383.0 24.33 465.9 577.7 897.7 559.9 527.7 477.5 -547.5 1662.0 334.5 468.5 406.2 342.0 1312.0 -117.1 79.507 1935.7 164.25 677.32 3150.86 1529.52 -923.282 24.432 269.67 407.47 1058.31 860.51 70.38 1260.02 283.8202 360.8016 656.528 1094.303 468.398 1859.374 4041.479
    2.0 -95.41801 0.0 174.1 117.3 2649.0 -628.07 778.3 1207.0 179.8 476.25 980.74 309.8 733.3 1857.0 224.8 165.3 2800.0 13.502 634.85 -2026.1 498.9 -44.958 -204.51 616.62 -56.69 417.6 NaN 46.06 NaN 470.4 NaN NaN -19.82 642.27 -174.6 179.7 967.9 141.1 388.4 NaN -339.8 2.406 -322.1 NaN 389.28 491.23 673.2271 249.18 1171.32 -86.231 9.389999 -456.08 775.56 324.316 -44.3408 566.2915 -92.7302 1230.781 1625.909 2979.973 -740.2421 2487.632 -2399.738
    3.0 16.07 -157.2 0.0 139.2 3972.0 1604.3 792.0 1356.0 146.2 -365.5 -274.54 170.5 -87.08 2036.0 139.67 -71.4 1044.7 -1553.9 -17.44 69.561 73.04601 133.66 66.214 269.0 595.2 134.1 936.3391 3736.0 77.083 331.6 -148.14 -82.28 -248.2 NaN 347.6 NaN 602.1 808.0 -81.68999 123.26 -126.2 134.6 -26.852 -1172.0 380.02 313.79 1083.63 92.5 -103.146 1412.0 392.56 849.08 -522.32 454.316 277.4415 825.9189 -1304.107 626.9603 161.7745 1690.065 -630.416 1014.686 4899.846
    4.0 47.2 -113.1 -45.33 0.0 3989.0 436.21 1050.2 1375.0 1001.0 683.6 -242.5 136.98 -595.1 2977.0 1250.0 -2631.0 4000.0 135.3 208.1 1352.5 -46.994 213.85 192.52 -106.2 -113.6 1358.0 391.044 2586.0 NaN 157.9 1856.45 69.0 277.0 NaN 88.93 NaN 234.2 -172.2 444.29 29.57 2303.0 -107.1 -26.486 -514.79 297.73 72.26 1696.96 76.73 -182.325 1000.8 164.05 558.08 -307.42 629.463 930.7075 -313.9384 -935.3928 524.3056 453.4621 1077.47 -478.068 1113.436 335.1104
    5.0 1606.0 1566.0 3049.0 2673.0 0.0 346.31 -801.9 83.91 -250.0 -281.4 973.8 235.9 816.7 -923.7 -355.1 104.6 -1114.0 -3061.2 123.5 -1295.0 238.1 -126.0 1314.8 925.6 1862.0 741.8 2100.0 NaN NaN 738.4 499.8 838.8 699.7 -148.9 190.4 1117.0 439.4 848.6 1036.0 NaN 403.8 3121.0 401.89 NaN -32.643 1201.3 186.24 1627.81 -150.961 -3745.0 495.28 NaN 34.94 1105.94 292.6468 1026.29 -348.7564 -2448.212 1178.048 -159.8631 -367.321 225.0235 -352.184
    6.0 82.593 -96.297 13.733 145.54 -1218.2 0.0 -328.5 -867.0 86.439 -392.5 299.23 220.7 -87.48 -495.25 -1508.5 -1039.0 -2012.0 -341.34 97.97301 -733.07 -16.521 -85.926 -139.58 -40.13 3000.0 374.2 NaN -332.4 75.71 -369.8 33.19 180.5 28.95 NaN -185.9 164.0 -43.88 -99.57999 101.2 NaN 308.7 68.972 -308.7 NaN -242.6 NaN -119.08 210.91 197.395 673.72 137.95 NaN 166.42 34.7884 -704.2285 -348.9494 -1230.497 -671.2604 53.7776 1741.061 1216.71 1037.56 NaN
    7.0 -17.253 -1301.0 332.3 24.144 1460.0 -524.3 0.0 -2686.0 190.5 -1545.0 -433.288 140.71 177.665 798.5 1524.0 274.5 158.4 -3178.5 -634.1 -1795.2 86.68999 134.1 NaN 1008.0 -1895.0 -595.7 NaN NaN NaN -123.8 372.5 NaN 822.2 NaN 117.0 419.8 NaN NaN -494.2 NaN 676.0 274.37 -75.7467 -804.28 509.3 659.22 934.67 64.43999 -1119.8 NaN 499.44 NaN -302.93 1038.3 659.2109 2770.851 NaN 2707.161 -1616.174 318.2128 4878.28 206.2641 NaN
    8.0 1987.0 191.6 2340.0 1825.0 465.4 265.5 148.4 0.0 -145.2 5.604 -212.9 NaN -329.3 NaN NaN NaN 542.0 -4080.9 NaN 401.88 NaN NaN NaN 2356.0 555.5 NaN NaN NaN NaN NaN -309.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 391.2 NaN NaN NaN NaN NaN NaN NaN 4911.4 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    9.0 199.0 91.811 -57.53 -146.6 653.3 394.78 770.6 -666.8 0.0 197.6 -16.486 -83.57 3645.0 NaN -47.97 -389.6 1732.0 -1908.7 -191.0 -109.51 -99.97601 -18.695 810.17 -208.71 1297.0 -35.89 -169.6 419.9 2.714 -986.0 478.5 346.6 -717.76 -62.43 -76.87 NaN 64.01 NaN 80.79201 NaN 64.21 437.739 -62.857 NaN -497.98 NaN 502.1 37.37 410.655 -153.7 NaN NaN NaN NaN 173.4451 824.1235 99.5504 147.2528 331.1518 571.6566 NaN 519.4855 NaN
    10.0 256.21 202.49 1011.0 1963.0 1590.0 -158.4 512.6 -410.21 -93.07999 0.0 -208.4 -160.7 209.0 NaN NaN NaN NaN NaN NaN 435.64 985.7 -41.537 NaN NaN NaN NaN NaN NaN -373.7 -742.7 NaN 114.3 NaN 236.6 NaN NaN -43.56 NaN 2371.0 NaN NaN 716.7 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    11.0 632.22 -582.82 622.73 1624.0 310.4 294.76 311.974 -224.4 33.415 389.7 0.0 342.4 195.3 NaN 3168.0 152.8 -1355.0 87.6 -193.23 62.031 -49.339 168.17 NaN -5.71 3351.0 9.222 NaN 861.1 NaN 80.68999 -72.07 82.96 -386.3 NaN 296.8 -92.12 -201.4 NaN 96.77 NaN -338.8 374.1 -28.231 745.4 -579.11 NaN 758.13 611.47 2973.58 NaN NaN NaN -1336.2 -1366.59 NaN NaN NaN NaN NaN NaN NaN NaN NaN
    12.0 238.5 -28.63 108.3 377.26 839.6 444.7 53.28 NaN 101.3 226.6 -251.7 0.0 NaN NaN NaN NaN NaN -9.2978 92.21 NaN NaN NaN NaN -142.2 1894.0 NaN NaN NaN 161.8 NaN NaN NaN NaN NaN NaN NaN 745.4 NaN 580.3 NaN NaN 227.6 NaN 489.15 260.64 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    13.0 -9.654 -844.3 179.0 375.0 650.9 475.2 433.207 -80.58 695.8 235.7 824.2 NaN 0.0 NaN NaN NaN NaN -199.94 1987.0 521.48 -208.6 492.9 -607.35 -425.4 974.0 -305.1 NaN 35.02 102.6 513.7 NaN -104.8 NaN -137.87 NaN NaN -422.7 155.7 NaN NaN NaN 397.0 -124.33 -454.92 -515.93 NaN NaN 1321.52 3858.07 NaN NaN NaN -209.209 NaN NaN NaN NaN NaN NaN NaN -584.728 NaN NaN
    14.0 326.04 498.8 -121.0 -45.44 -75.63 -467.95 -980.6 NaN NaN NaN NaN NaN NaN 0.0 1517.0 -472.4 NaN NaN -412.38 NaN NaN NaN NaN -65.76 2553.0 NaN NaN NaN -205.1 NaN NaN NaN NaN NaN NaN NaN NaN NaN 162.14 NaN NaN 124.3 -143.07 NaN NaN NaN 413.11 NaN NaN NaN NaN NaN 190.882 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    15.0 207.26 -124.32 105.63 -316.22 -660.2 -278.09 -851.0 NaN 119.5 NaN 3329.0 NaN NaN -1074.0 0.0 402.6 NaN NaN 242.2 NaN NaN NaN NaN -3.28 3888.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -330.2 -88.49 NaN NaN -131.9 -186.98 NaN NaN NaN 444.67 NaN NaN NaN NaN NaN 1178.79 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    16.0 205.65 -131.5 16.29 978.3 1876.0 39.33 -446.0 NaN 2831.0 NaN 160.8 NaN NaN 836.6 -639.9 0.0 NaN NaN NaN NaN NaN -473.0 NaN 215.9 1622.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -7.532001 965.0 NaN NaN 324.15 NaN NaN -420.24 NaN -729.5 2111.24 NaN NaN NaN NaN 2662.02 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    17.0 2257.3 3982.0 154.39 3969.0 1325.0 251.2 -131.0 -131.1 1460.0 NaN 3499.0 NaN NaN NaN NaN NaN 0.0 74.285 393.9 NaN 582.1 NaN NaN 3986.0 3709.2 NaN 3770.0 NaN NaN NaN 1268.0 NaN NaN NaN NaN NaN NaN NaN -391.9 NaN NaN 1371.0 NaN NaN NaN NaN 1271.42 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    18.0 -436.15 -13.317 1810.8 1698.1 -643.09 -230.38 -393.18 -41.594 307.16 NaN 45.309 -123.73 -430.49 NaN NaN NaN -558.39 0.0 746.12 -451.49 -67.106 -978.25 NaN 229.66 NaN NaN NaN NaN NaN NaN NaN NaN -713.65 NaN NaN NaN NaN 197.95 NaN NaN NaN 268.23 -250.25 NaN -516.31 NaN 489.93 NaN NaN -674.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    19.0 293.81 -181.93 111.8 170.1 -46.0 615.01 509.6 NaN 79.08 NaN 139.55 -33.64 -588.8 2412.2 -131.9 NaN 2987.0 -2307.8 0.0 NaN 176.5 -78.96 NaN 65.82 1283.0 117.53 NaN 468.8 -18.8 NaN 506.6 NaN -211.2 11.65 NaN 267.1 61.96 NaN 57.08 NaN -75.67 256.2 -28.653 NaN 237.42 NaN NaN NaN 297.646 NaN NaN NaN -50.8443 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    20.0 2017.7 -347.5 613.32 29.747 1525.8 1075.5 624.97 281.08 178.22 -188.0 59.594 NaN -310.82 NaN NaN NaN NaN -2617.7 NaN 0.0 27.618 94.606 NaN 701.95 -1398.7 NaN NaN NaN NaN NaN NaN 146.06 -18.328 NaN NaN NaN -447.95 NaN -421.21 NaN -271.176 1060.0 720.45 -65.631 508.72 NaN -369.31 -400.86 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    21.0 -65.68501 -359.6 -58.972 113.07 2177.0 1831.2 313.3 NaN 55.27 -888.3 48.852 NaN 872.0 NaN NaN NaN -338.0 -1592.8 -368.7 702.4 0.0 70.79 592.4 16.34 3985.0 24.44 1248.0 295.9 NaN 666.0 NaN NaN 128.8 NaN NaN NaN 280.0 NaN -70.45 NaN NaN -31.42 -325.77 530.3 207.12 NaN NaN NaN 955.8281 NaN NaN NaN 842.486 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    22.0 311.55 55.881 -142.2 -75.01 2389.0 1904.4 748.2 NaN -218.94 354.71 -461.35 NaN 215.3 NaN NaN 406.8 NaN 946.79 14.76 425.97 -66.21 0.0 187.43 46.29 3353.0 822.4 NaN NaN NaN -174.6 NaN 132.7 -139.6 NaN -178.3 NaN 160.7 NaN -147.64 NaN NaN 10.7 108.83 NaN 7.3664 NaN NaN 102.8 -1017.26 NaN NaN NaN 262.138 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    23.0 1302.6 255.41 -78.116 -38.939 963.37 893.38 NaN NaN -48.641 NaN NaN NaN 97.12801 NaN NaN NaN NaN NaN NaN NaN 603.29 1468.9 0.0 -323.17 NaN NaN NaN NaN NaN NaN NaN NaN 599.82 NaN NaN NaN 325.81 NaN NaN NaN NaN 289.08 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 401.698 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    24.0 -148.07 -663.45 -305.5 107.8 3139.0 2150.0 1282.0 2157.0 155.73 NaN 223.4 465.8 641.2 333.9 43.83 -825.9 2626.0 -2111.0 357.6 213.34 95.05 46.03 350.92 0.0 -131.8 441.5 3286.0 9.362 NaN 750.2 NaN 49.51 203.2 NaN 325.2 902.0 220.6 197.4 512.7 NaN NaN -37.183 190.45 NaN 22.779 NaN NaN NaN 516.966 1366.3 NaN NaN 120.59 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    25.0 3264.0 215.5 1885.0 -69.23 3664.0 2955.0 591.6 1554.0 1375.0 NaN -788.6 18.79 381.1 3873.0 -868.8 -94.87 1583.8 NaN 2331.0 1000.0 15.62 368.6 NaN 972.1 0.0 3986.0 -184.5 NaN NaN NaN NaN NaN -69.88 NaN NaN NaN NaN NaN -12.71 NaN NaN 207.16 96.855 NaN -27.161 NaN NaN NaN 593.203 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    26.0 -396.5 452.2 -330.0 1014.0 1091.0 1079.0 882.6 NaN -32.6 NaN -50.36 NaN 319.6 NaN NaN NaN NaN NaN -128.21 NaN 142.1 -423.1 NaN -65.74 3638.0 0.0 85.6 68.87 NaN NaN NaN 643.8 9.258 -70.24 NaN NaN 159.0 NaN 606.9 NaN NaN 119.3 53.75 NaN NaN NaN NaN NaN NaN 324.62 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    27.0 1744.0 NaN 1866.23 1931.39 316.6 NaN NaN NaN -328.1 NaN NaN NaN NaN NaN NaN NaN 1.655 NaN NaN NaN 1295.0 NaN NaN 167.5 2926.02 986.0 0.0 NaN NaN NaN NaN NaN 505.4 NaN NaN NaN NaN NaN NaN NaN NaN 2004.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    28.0 72.12 70.88 574.6 -1494.0 NaN 2645.0 NaN NaN 315.3 NaN 280.0 NaN 198.5 NaN NaN NaN NaN NaN 434.8 NaN -137.7 NaN NaN 52.01 NaN 655.7 NaN 0.0 NaN NaN NaN 212.4 NaN NaN NaN NaN -93.31 NaN NaN NaN NaN 29.45 166.56 NaN 89.744 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    29.0 -59.9 NaN -2.1662 NaN NaN 1334.0 NaN NaN 64.41 -397.5 NaN 13.97 -210.1 244.4 NaN NaN NaN NaN 41.54 NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN 467.1 NaN NaN NaN 356.6 NaN NaN -7.465 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    30.0 210.5 -146.1 0.4086 181.2 616.5 662.0 501.4 NaN 277.0 -214.8 -136.3 NaN -299.6 NaN NaN NaN NaN NaN NaN NaN -390.6 106.3 NaN 100.5 NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN -277.6 NaN NaN NaN NaN 96.59 NaN 778.78 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    31.0 28.17 NaN -93.18999 1041.33 -468.8 -3.428 -368.8 191.7 -72.58 NaN 69.25 NaN NaN NaN NaN NaN -818.8 NaN 11.72 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN -228.4 NaN NaN NaN 373.8 NaN NaN NaN NaN NaN NaN NaN -210.34 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    32.0 166.0 NaN 298.9 448.8 774.7 1965.0 NaN NaN -182.0 293.5 11.62 NaN 464.0 NaN NaN NaN NaN NaN NaN 780.71 NaN -23.81 NaN 186.4 NaN 17.81 NaN 200.6 NaN NaN NaN 0.0 -536.2 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -47.772 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    33.0 -62.08 160.4 187.5 -156.7 1439.0 648.8 -17.99 NaN 319.69 NaN 248.3 NaN NaN NaN NaN NaN NaN -763.19 362.6 753.21 -92.68 96.4 -364.76 -1360.0 981.5 121.4 81.44501 NaN NaN NaN NaN 558.0 0.0 NaN -83.7 NaN NaN NaN NaN NaN -378.1 -122.5 -186.4 NaN NaN NaN NaN NaN NaN NaN NaN NaN 2459.01 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    34.0 -22.04 -197.06 NaN NaN 1255.0 NaN NaN NaN -481.2 -93.066 NaN NaN 169.27 NaN NaN NaN NaN NaN -1428.0 NaN NaN NaN NaN NaN NaN 132.2 NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN 1025.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    35.0 477.1 154.0 -345.6 178.6 -452.3 145.0 -370.8 NaN 38.06 NaN -337.1 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -39.45 NaN -60.89 NaN NaN NaN NaN -360.0 NaN -47.81 NaN -116.7 NaN 0.0 NaN NaN NaN -133.35 NaN NaN NaN NaN NaN -322.46 NaN NaN NaN -3.44708 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    36.0 -291.9 -113.8 NaN NaN 1072.0 135.9 276.9 NaN NaN NaN 503.5 NaN NaN NaN NaN NaN NaN NaN -144.7 NaN NaN NaN NaN -194.9 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 -74.88 NaN -110.34 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    37.0 3.202 -300.6 1887.0 -115.9 959.3 2421.0 NaN NaN -153.4 945.6 -320.0 -479.1 -326.4 NaN NaN NaN NaN NaN -19.1 283.64 -207.3 -135.9 -199.87 -134.4 NaN 108.4 NaN 319.4 NaN 1168.0 NaN NaN NaN NaN NaN 1004.0 0.0 NaN 967.74 NaN -211.1 -24.82 321.62 NaN 185.82 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    38.0 -160.3 -139.5 -1367.0 -13.15 1253.0 1235.0 NaN NaN NaN NaN NaN NaN -528.8 NaN 904.1 35.16 NaN -1052.5 NaN NaN NaN NaN NaN -98.98 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 190.06 -57.38 NaN 112.7 -22.572 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    39.0 151.0 -152.2 58.76 64.28 -366.4 -269.7 -121.8 NaN 1955.5 -225.3 16.69 -285.5 NaN -112.76 -230.55 -311.9 650.7 NaN -160.48 93.773 -59.29 75.45 NaN -168.4 472.49 -340.9 NaN NaN -247.6 NaN -231.6 NaN NaN -416.5 221.74 97.36 -306.22 -229.97 0.0 NaN NaN 141.2 53.871 -310.13 NaN NaN NaN NaN -18.3768 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    40.0 -484.3 NaN 419.04 153.32 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 110.4 NaN 0.0 NaN 165.6 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    41.0 -314.6 698.5 670.8 -736.8 703.4 678.1 808.4 NaN -148.3 NaN 3.924 NaN NaN NaN NaN NaN NaN NaN 26.8 218.974 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 13.78 NaN NaN NaN 516.5 NaN NaN NaN 0.0 683.3 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    42.0 170.9 60.2 -2.619 191.5 2601.0 2540.7 1632.9 3630.0 364.423 1161.0 460.8 509.0 -214.1 621.9 1248.0 -198.32 1091.0 -5894.1 1336.0 578.3 370.6 224.4 -69.60101 60.78 223.782 522.9 2600.0 92.4 439.73 846.7 NaN NaN 476.9 NaN NaN NaN 81.56 21.04 666.5 -109.0 865.0 0.0 242.49 NaN 183.79 298.46 2187.67 499.19 -49.6851 313.43 NaN 373.49 -581.16 509.274 -883.9302 493.7747 -501.6695 -381.9636 -1078.495 1654.642 246.348 449.5142 NaN
    43.0 186.71 1182.6 47.23 199.48 -238.36 952.24 717.485 NaN 80.038 NaN 36.948 NaN 561.14 182.58 295.07 NaN NaN -1269.7 56.754 -140.77 70.075 -358.57 NaN -131.87 2991.9 -47.089 NaN 1.0902 NaN NaN NaN NaN 265.42 NaN NaN NaN 713.9 -7.56 -54.26 NaN NaN 20.834 0.0 NaN -523.96 NaN NaN NaN NaN NaN NaN NaN 19.221 NaN -617.0431 826.0889 NaN 580.1045 277.9974 -312.4461 NaN NaN NaN
    44.0 -21.23 NaN -1141.6 291.65 NaN NaN 594.45 NaN NaN NaN -447.04 -441.01 310.75 NaN NaN NaN NaN NaN NaN -14.016 17.052 NaN NaN NaN NaN NaN NaN NaN NaN -384.29 NaN NaN NaN NaN NaN NaN NaN NaN -367.48 NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    45.0 -44.069 -174.41 -201.52 -248.3 2985.8 4519.3 -523.8 NaN 945.14 NaN 966.35 -597.09 1368.0 NaN NaN -1035.8 NaN -1646.8 -642.44 -386.93 -175.29 -1.6641 NaN 14.947 4235.3 NaN NaN 40.987 NaN NaN NaN 92.429 NaN NaN 67.069 NaN -139.0 NaN NaN NaN NaN -61.922 1414.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN 184.792 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    46.0 -249.85 -734.87 -258.12 763.57 -516.99 NaN -588.21 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 499.59 NaN NaN NaN 0.0 NaN NaN NaN -1423.7 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    47.0 189.911 92.7 160.13 340.91 1460.82 -52.69 -649.31 NaN 141.01 NaN 32.71 NaN NaN -293.93 -196.23 -66.61 -2476.18 -2250.64 NaN -0.29 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 106.79 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 165.66 NaN NaN NaN NaN 0.0 -91.65 -188.913 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    48.0 82.6 -81.79 -14.89 246.08 -808.4 -132.93 -439.58 NaN -116.4 NaN -1031.78 NaN -361.25 NaN NaN -235.61 NaN NaN NaN -244.69 NaN -446.86 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 148.66 NaN NaN NaN NaN 117.57 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    49.0 183.02 90.0134 45.4296 119.767 868.681 -195.127 -128.903 NaN -201.105 NaN 81.6384 NaN -224.782 NaN NaN NaN NaN NaN -116.794 NaN -151.036 894.509 NaN -198.677 88.32201 NaN NaN NaN NaN NaN NaN NaN NaN NaN -72.16901 NaN NaN NaN -2.27584 NaN NaN 159.179 NaN NaN NaN NaN 308.235 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    52.0 -16.034 -39.373 -660.25 -139.78 7712.2 2676.7 NaN 2332.5 2311.5 NaN NaN NaN NaN NaN NaN NaN NaN 4998.6 NaN NaN NaN NaN NaN -580.46 NaN 11.442 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -30.564 NaN NaN NaN 1738.5 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    53.0 -41.0 57.86 -260.22 102.04 700.05 524.52 952.57 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN 912.2161 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    55.0 -0.067 298.55 -358.06 562.34 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 75.45 NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 -9.63211 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    56.0 438.76 -65.66 562.72 1661.26 942.22 1683.48 652.86 NaN NaN NaN 1915.2 NaN 1163.91 14.9127 1667.89 -1780.41 NaN NaN 180.813 NaN -313.202 -229.081 157.634 261.15 NaN NaN NaN NaN NaN NaN NaN NaN -271.226 NaN NaN NaN NaN NaN NaN NaN NaN 660.47 117.208 NaN -76.1382 NaN NaN NaN NaN NaN -604.897 125.526 0.0 200.04 NaN NaN NaN NaN NaN NaN NaN NaN NaN
    61.0 -309.943 -156.458 -388.19 -305.273 1343.95 641.041 -270.53 NaN NaN NaN 2622.04 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -108.172 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 742.31 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN
    84.0 77.7793 1753.135 -223.7499 365.3965 147.0663 152.7192 151.31 NaN 958.3726 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 360.0808 -323.6833 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 -103.0519 NaN -107.3735 NaN 133.3952 NaN -113.0145 4542.127
    85.0 996.5333 -1213.93 20.2719 -134.6708 -792.988 2127.895 1996.724 NaN 670.0023 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 986.9591 -450.7675 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 283.1673 0.0 -1403.31 NaN 855.7144 NaN NaN NaN NaN
    87.0 907.3842 -1151.296 26.3894 4101.374 1003.498 1142.897 NaN NaN -1504.191 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -935.6061 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -516.843 0.0 NaN NaN -1138.034 NaN NaN NaN
    89.0 454.5569 4896.219 2116.398 1259.47 -234.3916 985.4429 -399.1902 NaN 177.8232 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 2358.513 4927.307 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 1094.895 NaN NaN 0.0 -1679.977 NaN NaN NaN NaN
    90.0 -307.8627 159.8103 785.4488 3997.18 1044.876 94.629 275.7479 NaN -30.955 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 1039.753 122.9861 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -254.6989 NaN 4209.48 0.0 NaN NaN NaN NaN
    91.0 243.462 724.9113 -15.3556 -998.0735 662.1286 -558.6751 385.0068 NaN -1449.573 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -171.2349 186.7351 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -591.0368 NaN 1824.841 NaN NaN 0.0 NaN NaN NaN
    93.0 2660.16 232.905 -72.2694 15.0501 34.7132 -386.355 -444.615 NaN NaN NaN NaN NaN 178.501 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 1120.89 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN
    98.0 1131.305 1804.025 877.6001 -281.4692 448.381 1421.557 -371.2031 NaN 3306.365 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 4258.272 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -35.932 NaN NaN NaN NaN NaN NaN 0.0 NaN
    99.0 3504.325 658.1962 -969.6266 3792.062 1751.43 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -4116.333 NaN NaN NaN NaN NaN NaN NaN 0.0
]

Bnm = [
    0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 16.0 17.0 18.0 19.0 20.0 21.0 22.0 23.0 24.0 25.0 26.0 27.0 28.0 29.0 30.0 31.0 32.0 33.0 34.0 35.0 36.0 37.0 38.0 39.0 40.0 41.0 42.0 43.0 44.0 45.0 46.0 47.0 48.0 49.0 52.0 53.0 55.0 56.0 61.0 84.0 85.0 87.0 89.0 90.0 91.0 93.0 98.0 99.0
    1.0 0.0 -0.27232 0.0933 -0.4538 -4.674 -3.0099 -3.6156 -0.9977 0.1473 0.0 1.9294 -0.6215 -0.3155 4.9683 0.066729 1.857 -0.14836 -17.171 0.7335 -3.2647 -0.72772 1.2561 4.5311 -1.7109 15.89 -13.2 -2.693 1.521 -0.8557 0.9384 0.0 -0.3564 -0.499 -0.6581 3.457 0.0 0.009102 -1.042 0.6525 -1.679 -3.643 0.5481 0.70892 0.0 -0.26834 -2.0066 -3.1003 -6.2025 7.80959 0.2915 0.0 -0.2037 -1.1856 -1.96009 0.24511 -4.55446 -3.9547 2.05302 -4.2237 -1.11888 -2.55668 4.06308 -3.90698
    2.0 0.061708 0.0 -0.5886 -0.8552 -6.508 10.0 0.1482 -1.955 0.69911 0.0 -2.4224 0.0 -2.509 -8.653001 0.0 0.0 -10.72 0.0 0.0 8.1549 -1.487 0.0 0.0 -2.1164 9.805 0.08726 NaN 1.545 NaN 0.0 NaN NaN 0.5941 -3.8641 1.96 0.0 -0.03862 -0.3025 0.0 NaN 1.297 -0.1882 -0.20367 NaN 0.0 -1.8285 1.9523 0.0 3.92127 0.0 0.0 2.4334 0.0 0.0 -0.65212 -0.14057 -0.25088 0.17033 -2.40107 4.82369 0.0 -3.28182 -0.27451
    3.0 -0.2998 0.6166 0.0 -0.65 -13.16 -2.0299 -1.726 -2.118 -1.237 1.874 0.91491 -0.02393 -0.1859 -8.729 0.037692 0.7078 -1.7112 0.1615 0.9437 1.8881 -0.21322 -0.46141 -0.63629 -1.776 3.309 -0.8156 -4.4141 -25.0 -0.3094 0.3778 2.2658 0.5677 1.214 NaN -1.43 NaN -7.798 -5.331 1.1229 -0.6133 0.386 -1.231 -0.44215 10.106 -0.23326 -1.1552 0.9144 1.4186 2.6565 -3.9917 0.0 -6.5044 3.8856 -1.21943 -0.1843 -4.35541 3.87435 -4.30606 -0.56713 -3.38808 0.0 -1.73087 -3.53433
    4.0 0.3575 1.172 0.4223 0.0 -14.09 1.9094 -1.9939 -1.702 -1.871 -1.02 2.72 0.0 2.978 -19.16 0.0 13.56 -16.68 -2.369 1.333 0.0 0.30443 0.19998 -0.25363 0.6081 19.72 -9.968 0.196357 -9.5 NaN 1.408 5.9959 0.4317 -0.5955 NaN 0.2236 NaN -0.2311 0.2257 0.1947 0.0 -6.346 0.2564 -0.29519 4.9372 -0.24204 -0.1919 -0.7064 0.9394 3.36729 -2.5194 0.0 0.565 4.1819 -1.39039 -4.30401 2.06128 1.26386 -4.26859 -2.03845 -2.43883 0.0 -0.428312 2.26653
    5.0 -4.746 -5.809 -12.77 -5.765 0.0 -2.4583 3.824 -1.262 2.857 2.379 -5.633 -0.5874 -5.092 2.468 0.58 -5.014 5.916 9.5413 0.8503 4.3634 -0.7077 0.1322 0.0 -2.027 32.07 -2.098 0.0 NaN NaN -1.771 -2.41 0.0 -1.767 1.034 -1.2 -2.196 0.0 0.0 -2.995 NaN -0.9346 -13.69 -0.43627 NaN -0.010434 -4.6803 -3.1349 -5.8598 0.527951 4.5137 -2.6496 NaN 1.0075 -1.68763 -1.59145 -3.17777 1.92474 -2.13997 5.33766 2.41804 0.0 -4.3271 -1.58363
    6.0 -0.48575 0.6304 -0.11768 -0.48799 9.7928 0.0 1.0823 -1.258 -0.46505 2.256 -1.2702 -0.6402 -0.5522 1.0807 4.4917 5.603 13.46 -5.8042 -0.28666 2.3351 -0.28137 -0.26371 0.0 -0.06709 -11.81 -2.212 NaN 4.419 -0.3753 2.203 -0.2074 0.3161 -0.3303 NaN 0.03906 0.0 0.0 0.0 -1.078 NaN -1.175 -0.41997 1.7454 NaN 0.29562 NaN 0.9938 -2.4953 -0.985958 -2.3567 -1.4112 NaN 0.994 0.0 -1.26236 1.41462 3.91001 1.58636 -0.9791 -4.73718 -5.44835 -5.46099 NaN
    7.0 0.8389 4.072 1.158 1.6504 -8.673 4.6065 0.0 19.44 -3.669 6.512 3.0862 0.0 -3.72906 -5.869 -2.531 -0.5905 0.5246 9.8213 3.59 12.708 1.592 0.9495 NaN -1.795 9.303 2.634 NaN NaN NaN 3.847 -0.9091 NaN 0.0 NaN -0.611 2.436 NaN NaN 2.826 NaN -0.9909 -0.5861 -0.98514 2.828 -0.29 -1.8841 -4.4232 -0.0094 7.95909 NaN -6.6553 NaN 0.5975 0.0 -0.59252 -2.00944 NaN -2.51124 5.32305 -1.42376 -6.38758 -2.57259 NaN
    8.0 -4.615 0.4936 -5.043 -3.743 -1.841 -2.905 -2.757 0.0 -0.738 0.0 0.0 NaN 0.0 NaN NaN NaN -6.792 12.859 NaN 0.0 NaN NaN NaN -3.347 0.0 NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -0.8799 NaN NaN NaN NaN NaN NaN NaN 9.8501 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    9.0 -0.8709 -0.71715 1.212 0.2419 -1.412 -0.36048 -0.5873 1.918 0.0 0.0 -0.27924 0.0 -26.91 NaN 0.0 0.1944 -9.896001 4.1057 0.6835 0.96888 0.0 -0.52606 -3.2209 0.8293 -1.365 -0.1009 2.986 0.09772 0.0 3.847 -1.148 -0.7017 2.9431 -0.8977 0.5372 NaN -1.077 NaN -3.702 NaN 0.0 -2.7983 0.2898 NaN 0.79715 NaN 0.3401 0.0 -1.1056 0.242 NaN NaN NaN NaN -5.52695 -2.88368 0.45575 -2.74448 -1.11421 -4.43975 NaN -5.4673 NaN
    10.0 0.0 0.0 -2.167 2.656 -24.57 -0.6469 -2.145 0.0 0.0 0.0 0.0 0.0 -0.6241 NaN NaN NaN NaN NaN NaN 0.0 -2.986 0.0 NaN NaN NaN NaN NaN NaN 0.0 0.0 NaN 0.0 NaN 0.0 NaN NaN 0.0 NaN 0.0 NaN NaN -1.516 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    11.0 -3.3912 1.6732 -1.7605 -9.409 1.538 0.3745 -1.3412 0.0 0.21913 0.0 0.0 0.0 -9.75 NaN -24.07 -1.099 7.683001 1.4109 0.43006 1.0567 0.0 -1.0536 NaN -0.2724 -14.54 -0.3292 NaN -1.693 NaN 0.0 0.0 0.0 1.892 NaN -1.264 0.9031 0.5487 NaN 0.4827 NaN 1.837 -1.976 0.0 0.0 0.9455 NaN -0.2593 -2.4783 -5.45555 NaN NaN NaN 6.7742 5.90673 NaN NaN NaN NaN NaN NaN NaN NaN NaN
    12.0 -0.5358 0.0 -0.262 0.0 -1.215 -0.06819 0.0 NaN 0.0 0.0 0.0 0.0 NaN NaN NaN NaN NaN 0.0 0.0 NaN NaN NaN NaN 0.572 0.0 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN -0.622 NaN 0.0 NaN NaN -0.7762 NaN 0.0 -1.2868 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    13.0 -0.03242 2.945 0.05615 -1.57 -0.7132 0.1198 -0.605276 0.0 -0.9619 0.1314 -6.009 NaN 0.0 NaN NaN NaN NaN -6.9867 -8.022 0.0 -0.2571 -6.475 2.3467 0.9514 -1.368 0.7063 NaN 1.702 0.0 0.0 NaN 0.0 NaN 0.0 NaN NaN 2.23 -1.578 NaN NaN NaN -1.379 -0.29398 0.0 0.3835 NaN NaN 0.0 -9.278001 NaN NaN NaN -0.470423 NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN
    14.0 -2.6348 -5.148 -1.901 5.141 -0.1511 0.65743 3.671 NaN NaN NaN NaN NaN NaN 0.0 -12.72 -0.2051 NaN NaN -0.49092 NaN NaN NaN NaN -0.3148 -11.9 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN -1.147 0.0 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    15.0 -1.0916 0.0 -0.60667 0.0 1.743 -0.39888 1.034 NaN 0.0 NaN -13.78 NaN NaN 9.0 0.0 -1.614 NaN NaN 0.0 NaN NaN NaN NaN 0.0 -16.26 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 0.0 NaN NaN -0.007369 0.0 NaN NaN NaN 0.8351 NaN NaN NaN NaN NaN -9.99911 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    16.0 -1.4436 0.0 -0.6022 -6.481 11.5 -3.354 -0.7738 NaN -13.01 NaN 0.8719 NaN NaN -5.208 2.561 0.0 NaN NaN NaN NaN NaN 0.8883 NaN -1.399 -4.812 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -0.8077 0.0 NaN NaN -1.4264 NaN NaN 0.26315 NaN 6.2811 0.0 NaN NaN NaN NaN -13.9855 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    17.0 -5.6676 -19.72 1.2458 8.497001 -6.263 -1.274 -0.7957 -0.9399 -13.73 NaN -22.96 NaN NaN NaN NaN NaN 0.0 14.613 -4.703 NaN -2.33 NaN NaN -16.15 0.0 NaN 1.686 NaN NaN NaN -3.042 NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN -3.323 NaN NaN NaN NaN -6.4373 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    18.0 3.4225 0.0 -5.7594 -9.8887 2.9346 -0.3825 0.5546 -3.5887 -0.8161 NaN -1.4281 0.0 -0.728 NaN NaN NaN 1.2436 0.0 -2.1937 -0.4031 -0.4933 2.1884 NaN -1.5005 NaN NaN NaN NaN NaN NaN NaN NaN 1.9191 NaN NaN NaN NaN -0.3316 NaN NaN NaN -1.4137 0.3811 NaN 1.6943 NaN -0.5591 NaN NaN 1.8854 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    19.0 -1.3979 0.0 -0.5959 -0.8218 4.415 -0.9444399 -1.916 NaN -0.3808 NaN -0.43665 0.0 1.481 -10.495 0.0 NaN -9.336 5.8186 0.0 NaN -1.237 0.0 NaN -0.6265 3.361 -0.44692 NaN 0.0 -0.3652 NaN -1.263 NaN 0.0 -1.229 NaN 0.0 -0.4161 NaN -0.7417 NaN 0.0 -1.526 -0.48152 NaN -1.2928 NaN NaN NaN -0.66206 NaN NaN NaN -0.025287 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    20.0 -9.0933 1.216 -1.595 0.0 -4.9155 -3.4339 -4.6878 0.0 -0.91676 0.0 -0.71198 NaN 0.0 NaN NaN NaN NaN 4.614 NaN 0.0 0.0 0.0 NaN -1.7576 0.0 NaN NaN NaN NaN NaN NaN 0.0 0.0 NaN NaN NaN 1.5141 NaN 0.0 NaN 1.4107 -2.822 -1.5187 0.0 -1.4005 NaN 0.0 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    21.0 0.074091 1.238 0.1046 -0.59569 -4.363 -2.9694 -0.5041 NaN 0.0 3.426 0.0 NaN -2.939 NaN NaN NaN 3.372 2.8397 1.992 0.0 0.0 0.0 -4.2459 0.07287 -15.7 -0.4713 -2.04 0.0 NaN 0.0 NaN NaN -0.2077 NaN NaN NaN 0.0 NaN 0.0 NaN NaN -0.248 2.0412 0.0 -0.43964 NaN NaN NaN -1.00012 NaN NaN NaN -1.8434 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    22.0 -1.1856 0.0 0.39662 -0.42684 -4.509 -2.7981 -1.319 NaN 1.0749 0.0 1.8569 NaN -2.482 NaN NaN -1.524 NaN -4.9886 0.0 0.0 0.0 0.0 -3.446 -0.2115 -14.2 -2.05 NaN NaN NaN 0.0 NaN -0.1183 1.022 NaN 0.7426 NaN 0.0 NaN 0.0 NaN NaN -0.2702 -0.86059 NaN 0.40456 NaN NaN -0.0879 3.0909 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    23.0 -8.427 0.0 0.69975 -0.19378 0.0 0.0 NaN NaN -0.79498 NaN NaN NaN -0.64386 NaN NaN NaN NaN NaN NaN NaN -3.977 -5.0 0.0 1.1973 NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN 0.0 NaN NaN NaN NaN -0.5852 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    24.0 1.0927 2.3281 2.12 -0.6785 -5.964 -2.813 -2.815 -3.718 1.4412 NaN 0.1237 -0.7730001 -1.486 -0.1415 0.0 3.15 -10.59 6.3625 0.7676 2.1861 -0.2348 0.05388 -1.3456 0.0 9.802 -0.5353 0.0 1.033 NaN 0.8165 NaN 0.6829 -1.328 NaN 0.4405 0.0 0.3756 -0.4858 0.02702 NaN NaN -0.047827 0.012716 NaN 0.042145 NaN NaN NaN 0.19471 -3.8385 NaN NaN 1.2816 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    25.0 -20.84 -1.519 -10.98 -0.7359 34.13 -6.383 -3.08 0.0 2.499 NaN 1.693 0.0 -5.682 -9.315001 2.948 -9.612 0.0 NaN -9.238001 0.0 -1.099 -1.748 NaN -6.82 0.0 25.0 2.2806 NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN -0.4758 NaN NaN -1.2206 -1.2993 NaN 0.611 NaN NaN NaN -0.7904 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    26.0 3.092 -1.997 3.066 -4.702 -1.274 -1.766 -2.606 NaN 0.0847 NaN 0.02448 NaN -1.368 NaN NaN NaN NaN NaN 0.50349 NaN -0.153 0.8154 NaN 0.0967 -14.25 0.0 0.0 3.217 NaN NaN NaN -0.7376 -0.1079 0.0 NaN NaN 0.0 NaN 0.0 NaN NaN -0.5691 0.0 NaN NaN NaN NaN NaN NaN -1.1828 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    27.0 -4.082 NaN -9.28137 -5.15065 0.0 NaN NaN NaN 0.3045 NaN NaN NaN NaN NaN NaN NaN -1.164 NaN NaN NaN -4.224 NaN NaN 0.0 -1.8141 0.0 0.0 NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN -4.5 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    28.0 -1.126 -1.09 -3.702 7.679 NaN -9.172001 NaN NaN -0.5617 NaN -0.9491 NaN -1.434 NaN NaN NaN NaN NaN 0.0 NaN 0.0 NaN NaN -0.9095 NaN -2.203 NaN 0.0 NaN NaN NaN -0.3692 NaN NaN NaN NaN -0.1286 NaN NaN NaN NaN -0.9194 -1.0407 NaN -1.0122 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    29.0 -0.008313 NaN 0.1196 NaN NaN -1.863 NaN NaN 0.0 0.0 NaN 0.0 0.0 0.0 NaN NaN NaN NaN 0.646 NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN 0.0 NaN NaN NaN 0.0 NaN NaN -0.4151 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    30.0 -1.081 0.0 -0.4601 -1.081 -0.1457 -1.903 -1.939 NaN -1.357 0.0 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN 0.0 0.0 NaN -0.8269 NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN -0.7691 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    31.0 0.0 NaN 1.5036 -5.2068 2.421 0.483 0.7775 0.0 0.4909 NaN 0.0 NaN NaN NaN NaN NaN 3.229 NaN 0.7004 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN 0.0 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN 0.9103 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    32.0 -0.7116 NaN -0.8374 -0.9919 0.0 -3.858 NaN NaN 0.7905 0.0 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN 0.0 NaN -0.09203999 NaN -0.7294 NaN -0.2245 NaN -0.228 NaN NaN NaN 0.0 1.995 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -0.009201 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    33.0 -0.3658 -0.5148 -0.9020001 0.1221 -1.673 0.8050001 0.0 NaN -1.4362 NaN -1.198 NaN NaN NaN NaN NaN NaN 0.603 0.0 0.0 -0.01307 -0.776 0.0 7.402 0.0 0.05397 0.0 NaN NaN NaN NaN -1.924 0.0 NaN 0.1436 NaN NaN NaN NaN NaN 2.56 0.0913 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    34.0 -0.1018 1.8723 NaN NaN -2.538 NaN NaN NaN 9.604 0.0 NaN NaN 0.0 NaN NaN NaN NaN NaN 7.698 NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    35.0 -2.175 -1.303 1.545 -0.7113 1.956 -0.6343 0.1043 NaN -0.5189 NaN 0.8843 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -0.8556 NaN -0.6321 NaN NaN NaN NaN 0.0 NaN 0.0 NaN -0.06775 NaN 0.0 NaN NaN NaN 0.1022 NaN NaN NaN NaN NaN -0.020828 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    36.0 0.0 0.0 NaN NaN -1.248 0.0 -0.9948 NaN NaN NaN -1.792 NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    37.0 -0.6894 -0.09576 -8.707001 -0.02925 0.0 0.0 NaN NaN 3.143 0.0 1.052 0.8031 0.24 NaN NaN NaN NaN NaN 1.195 1.5491 0.0 0.0 0.0 -0.3226 NaN 0.0 NaN 0.03419 NaN 0.0 NaN NaN NaN NaN NaN 0.0 0.0 NaN 0.0 NaN 0.0 -0.5324 -4.9963 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    38.0 -0.0108 0.9076 7.883 0.001763 0.0 0.0 NaN NaN NaN NaN NaN NaN 2.822 NaN 0.0 0.09722999 NaN 1.6854 NaN NaN NaN NaN NaN -0.02128 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 0.0 0.0 NaN -0.4103 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    39.0 -0.9023 0.0 -0.6263 -0.5839 1.129 1.732 -0.6029 NaN -10.0 0.0 -0.5966 0.0 NaN 0.0 0.0 0.0 0.0 NaN 1.0644 0.0 0.0 0.0 NaN -0.123 -0.5503 0.0 NaN NaN 0.0 NaN 0.0 NaN NaN 0.0 -0.5021 0.0 0.0 0.0 0.0 NaN NaN -0.8783 -0.67747 0.0 NaN NaN NaN NaN -0.450261 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    40.0 2.467 NaN -0.1703 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    41.0 1.287 -2.159 -1.807 3.332 -1.383 -1.335 -2.929 NaN 0.0 NaN -1.174 NaN NaN NaN NaN NaN NaN NaN 0.0 -0.8061 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -0.936 NaN NaN NaN 0.0 NaN NaN NaN 0.0 6.321 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    42.0 -0.8062 0.1565 1.094 -0.5561 -1.25 -3.5236 -2.8719 -6.23 2.1022 -0.5724 -0.06206 -0.2248 1.134 0.258 -2.79 1.452 -0.2657 13.088 -3.395 1.493 -0.322 -0.1399 -0.4274 0.024262 -0.5826 -0.2485 -4.8 1.191 -0.2619 0.2545 NaN NaN -0.61 NaN NaN NaN 0.7293 -0.1288 -0.1555 0.0 -2.89 0.0 -0.038323 NaN -0.15182 -0.6823 -1.057 0.5131 3.86633 0.0911 NaN 0.0647 4.1194 -0.533745 4.0287 -4.10779 1.6261 -3.00577 3.65154 -4.24837 -1.9006 2.1727 NaN
    43.0 -1.3546 -5.0 0.64039 0.033695 5.0 -3.3287 -1.48515 NaN -0.10124 NaN 0.0 NaN -0.7058 0.0 0.0 NaN NaN 2.6335 0.89781 0.309 -1.149 1.3307 NaN -0.014204 -9.3959 0.0 NaN 1.5927 NaN NaN NaN NaN 0.0 NaN NaN NaN -2.7759 0.0 1.0612 NaN NaN -0.34718 0.0 NaN 0.4945 NaN NaN NaN NaN NaN NaN NaN 0.4413 NaN 2.40181 -0.93908 NaN -2.24427 0.25973 1.126 NaN NaN NaN
    44.0 0.0 NaN 8.656199 -1.2039 NaN NaN -2.2535 NaN NaN NaN 0.0 0.0 0.0 NaN NaN NaN NaN NaN NaN 0.0 0.0 NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    45.0 0.027778 0.0 -0.068774 0.25471 -6.227 -19.456 3.158 NaN -3.8168 NaN -2.1861 2.5295 -2.6254 NaN NaN 3.078 NaN 2.7237 2.8574 2.3961 0.32745 -0.37825 NaN -0.11086 -16.954 NaN NaN 1.1526 NaN NaN NaN 0.17437 NaN NaN 0.23965 NaN 0.0 NaN NaN NaN NaN -0.059442 -2.8776 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    46.0 1.7054 3.3463 1.4084 -1.3961 2.1888 NaN 0.9707 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -0.8158 NaN NaN NaN 0.0 NaN NaN NaN 3.9162 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    47.0 -0.9413 -0.8135 -1.3578 -1.008 -5.962 -0.7465 1.3487 NaN -0.9829 NaN -0.8641 NaN NaN 0.0 -0.0511 -0.696 10.9894 5.4638 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -0.9536999 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -0.9151 NaN NaN NaN NaN 0.0 -0.3411 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    48.0 -0.615 0.0 -0.8126 -0.8198 3.8321 0.6584 0.3142 NaN 0.0 NaN 3.9906 NaN 0.0 NaN NaN 0.0 NaN NaN NaN 0.0 NaN 0.4276 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -1.0328 NaN NaN NaN NaN 0.7114 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    49.0 -1.22588 -0.920017 -0.570848 -0.608934 -2.20103 1.17652 -0.468519 NaN 0.7123 NaN -1.3111 NaN -0.410609 NaN NaN NaN NaN NaN 0.54754 NaN -0.458824 -2.93164 NaN 0.435084 0.1605 NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN 0.601128 NaN NaN -1.11479 NaN NaN NaN NaN 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    52.0 -0.5337 0.0 1.4602 0.1146 -7.6896 -5.4899 NaN -4.0832 -6.6054 NaN NaN NaN NaN NaN NaN NaN NaN 1.1962 NaN NaN NaN NaN NaN 1.1653 NaN 0.7866 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -0.2455 NaN NaN NaN -3.9729 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    53.0 0.0 0.0 0.0 0.0 0.2237 0.4314 -0.8437 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN -2.02314 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    55.0 -0.3935 -1.4153 3.6831 -3.4373 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -0.6572 NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 0.555403 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    56.0 -1.2256 0.0 -2.4316 -8.697001 -2.3183 -8.003099 -1.8013 NaN NaN NaN -9.9994 NaN -7.1332 0.0 -7.58704 5.71975 NaN NaN -0.0641231 NaN 0.4139 0.0 0.0 -0.8815 NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN -2.0132 -0.3946 NaN 0.0 NaN NaN NaN NaN NaN 1.04748 0.289226 0.0 0.1979 NaN NaN NaN NaN NaN NaN NaN NaN NaN
    61.0 0.64856 0.0 1.04305 0.669085 -2.16488 0.0 0.0 NaN NaN NaN -11.7013 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -0.132539 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -4.1234 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN
    84.0 0.83298 -3.112 0.2781 -1.20387 2.08434 2.17856 -3.00895 NaN -4.73012 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -1.37211 -0.13058 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 -1.42845 NaN -3.90856 NaN -2.64806 NaN -1.48378 0.85207
    85.0 -5.48462 5.14012 -0.33187 1.38462 3.87681 -4.14737 -5.47448 NaN -5.05233 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -2.85003 0.94786 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -5.14586 0.0 -4.02583 NaN -1.14021 NaN NaN NaN NaN
    87.0 0.46568 4.7125 0.38631 -1.19824 -2.48421 -2.92427 NaN NaN 4.63058 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 4.67403 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -3.5941 0.0 NaN NaN 0.70013 NaN NaN NaN
    89.0 -1.47382 -1.86077 -2.11062 -1.91791 4.01945 -2.42145 0.16195 NaN -2.24998 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -1.30922 -0.63511 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -4.10982 NaN NaN 0.0 4.40837 NaN NaN NaN NaN
    90.0 3.60996 1.85317 -3.09114 -2.994 -3.71273 0.30192 -1.97057 NaN -0.71244 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -3.51543 -1.96464 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -0.46519 NaN -0.75036 0.0 NaN NaN NaN NaN
    91.0 -1.03559 -4.79411 -0.15571 4.48792 -3.28127 1.11398 -2.49001 NaN 4.44554 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.62135 -0.9786 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.45605 NaN -4.26948 NaN NaN 0.0 NaN NaN NaN
    93.0 -6.86512 0.0 0.0 0.0 0.0 3.77874 -0.776483 NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -3.1028 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN
    98.0 -0.93175 -5.18107 -1.96581 1.00738 -1.61662 0.34143 0.99046 NaN 0.15196 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.34026 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -1.52136 NaN NaN NaN NaN NaN NaN 0.0 NaN
    99.0 2.66288 0.57714 3.81305 -2.55407 -1.09998 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -1.71395 NaN NaN NaN NaN NaN NaN NaN 0.0
]

Cnm = [
    0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0 8.0 9.0 10.0 11.0 12.0 13.0 14.0 15.0 16.0 17.0 18.0 19.0 20.0 21.0 22.0 23.0 24.0 25.0 26.0 27.0 28.0 29.0 30.0 31.0 32.0 33.0 34.0 35.0 36.0 37.0 38.0 39.0 40.0 41.0 42.0 43.0 44.0 45.0 46.0 47.0 48.0 49.0 52.0 53.0 55.0 56.0 61.0 84.0 85.0 87.0 89.0 90.0 91.0 93.0 98.0 99.0
    1.0 0.0 0.0 0.0 0.0 0.001551 0.0 0.001144 0.0 0.0 0.0 -0.0031331 0.0 0.0 -0.010252 0.0 0.0 0.0 0.036 0.0 0.009198 0.0 0.0 -0.008735 0.003388 -0.04831 0.02156 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.00098 -0.0020983 0.0 0.0 0.0 0.0 0.00975 -0.008808 0.0 0.0 -0.00016 0.0 0.0 0.0 0.0067607 0.0103753 -0.0166992 0.0074475 0.0 0.0 0.0 -0.0254482
    2.0 0.0 0.0 0.0 0.0 0.004822 -0.014972 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.01088 0.0 0.0 0.01339 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 NaN 0.0 NaN 0.0 NaN NaN 0.0 0.0055 0.0 0.0 0.0 0.0 0.0 NaN 0.0 0.0 0.004517 NaN 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0141503 0.0 -0.0272505 0.0 0.0 0.0329518
    3.0 0.0 0.0 0.0 0.0 0.01208 0.0 0.0 0.0 0.004237 0.0 0.0 0.0 0.0 0.008138 0.0 0.0 0.0 0.0082 0.0 0.0 0.0 0.0 0.0 0.002645 -0.02844 0.00145 0.012232 0.04593 0.0 0.0 -0.00353 0.0 0.0 NaN 0.0 NaN 0.01966 0.006077 -0.00105 0.0 0.0 0.001488 0.0 -0.01428 0.0 0.0 -0.002141 1.8e-5 -0.002819 0.0 0.0 0.01219 -0.004212 0.0 0.0 0.0085572 0.0039496 0.0 0.0 0.0202432 0.0 0.0 -0.0358313
    4.0 0.0 0.0 0.0 0.0 0.0153 0.0 0.0 0.0 0.000239 0.000869 -0.003449 0.0 0.0 0.03333 0.0 -0.007036 0.02112 0.0039 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.02093 2.24992e-5 0.008819 NaN 0.0 -0.01862 0.0 0.0 NaN 0.0 NaN 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.00045 -0.004625 0.0 0.0 -0.00331 -0.007723 0.0 0.007841 0.0 0.0041351 0.0 0.0 0.0 0.0 0.0 0.0239277
    5.0 0.0009181 0.005197 0.01435 -0.000332 0.0 0.0029287 -0.007514 0.0 -0.006022 -0.006668 0.00769 0.0 0.006065 0.0 0.0 0.008854 -0.007126 0.1842 -0.002478 0.0 0.0 0.0 0.0 0.0 -0.009397 0.0 0.0 NaN NaN 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 NaN 0.0 0.01446 -0.002004 NaN 0.0 0.0 0.007737 0.005497 0.0 0.0636 0.00159 NaN -0.002547 0.0 0.0 0.0 0.0 0.0332832 0.0 0.0 0.0 0.0030855 0.0090199
    6.0 0.0 -0.0018 0.0 0.0 -0.016158 0.0 -0.0022 0.02998 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.006551 -0.02004 0.0141 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.002688 NaN -0.007797 0.0 0.0 0.0 0.0 0.0 NaN 0.0 0.0 0.0 0.0 0.0 NaN 0.0 0.0 -0.00335 NaN -9.0055e-5 NaN 0.0 0.0048 0.0 0.0 0.00172 NaN -0.003346 0.0 0.0149627 0.0 0.0 -0.0048637 0.0 0.0 0.0 0.0 NaN
    7.0 0.0009021 0.0 0.0 0.0 0.01641 -0.004 0.0 -0.02702 0.008838 0.0 -0.002012 0.0 0.010763 0.01032 0.0 0.002205 0.0 -0.0034 0.0 -0.015455 0.0 0.0 NaN 0.0 0.0 0.0 NaN NaN NaN 0.0 0.0 NaN 0.0 NaN 0.0 0.0 NaN NaN 0.0 NaN 0.0 -0.00030011 0.0033318 0.0 0.0 0.0 0.006054 0.0 -0.008555 NaN 0.01656 NaN 0.0 0.0 0.0451602 0.0 NaN -0.0022558 0.0 0.0 0.0 0.0032526 NaN
    8.0 0.0 0.0 0.0 0.0 0.0 0.002283 0.002329 0.0 0.0 0.0 0.0 NaN 0.0 NaN NaN NaN 0.01655 -0.0101 NaN 0.0 NaN NaN NaN 0.0 0.0 NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    9.0 0.0 0.0 -0.003715 0.0001133 0.000954 0.0 -0.003252 0.0 0.0 0.0 0.0 0.0 0.04757 NaN 0.0 0.001863 0.0141 0.0 0.0 0.0 0.0 0.0 0.0021443 -0.0014 -0.02253 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 NaN 0.0 NaN 0.011586 NaN 0.0 0.00364 0.0 NaN 0.00029817 NaN 0.0 0.0 0.0 0.0 NaN NaN NaN NaN 0.0160979 0.0041569 0.0 0.0 0.0 0.0 NaN 0.0068391 NaN
    10.0 0.0 0.0 0.0 -0.01355 0.06212 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 NaN NaN NaN NaN NaN NaN 0.0 0.0 0.0 NaN NaN NaN NaN NaN NaN 0.0 0.0 NaN 0.0 NaN 0.0 NaN NaN 0.0 NaN 0.0 NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    11.0 0.0039282 0.0 0.0 0.01338 -0.004885 0.0 0.001074 0.0 0.0 0.0 0.0 0.0 0.04051 NaN 0.04303 0.0 -0.01012 0.0 0.0 0.0 0.0 0.0 NaN 0.0 0.0 0.0 NaN 0.0 NaN 0.0 0.0 0.0 0.0 NaN 0.0 0.0 0.0 NaN 0.0 NaN 0.0 0.001682 0.0 0.0 0.0 NaN 0.0 0.0 0.0 NaN NaN NaN -0.0082 -0.007325 NaN NaN NaN NaN NaN NaN NaN NaN NaN
    12.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 NaN 0.0 0.0 0.0 0.0 NaN NaN NaN NaN NaN 0.0 0.0 NaN NaN NaN NaN 0.0 0.0 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN 0.0 NaN 0.0 NaN NaN 0.0 NaN 0.0 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    13.0 0.0 0.0 0.0 0.0 0.000815 0.0 -0.000914 0.0 -0.002462 0.0 0.008271 NaN 0.0 NaN NaN NaN NaN 0.0193 0.01065 0.0 0.002418 0.01806 0.0 0.0 -0.01983 0.0 NaN 0.0 0.0 0.0 NaN 0.0 NaN 0.0 NaN NaN 0.0 0.0 NaN NaN NaN 0.0 0.0 0.0 0.0 NaN NaN 0.0 7.0e-6 NaN NaN NaN 0.006918 NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN
    14.0 0.0033576 0.01039 0.006999 -0.0142 0.0 0.0 -0.005908 NaN NaN NaN NaN NaN NaN 0.0 0.02557 0.01058 NaN NaN 0.0062553 NaN NaN NaN NaN 0.0 0.0 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN 0.0 0.0 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    15.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 NaN 0.0 NaN 0.01193 NaN NaN -0.01795 0.0 0.0 NaN NaN 0.0 NaN NaN NaN NaN 0.0 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 0.0 NaN NaN 0.0 0.0 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN 0.020758 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    16.0 0.0 0.0 0.0 0.007088 0.09 0.006714 0.002634 NaN 0.01558 NaN 0.0 NaN NaN 0.004801 0.0 0.0 NaN NaN NaN NaN NaN 0.0 NaN 0.0 -0.01856 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 0.0 NaN NaN 0.0 NaN NaN 0.0 NaN 0.0 0.0 NaN NaN NaN NaN 0.017738 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    17.0 0.0 0.02783 0.0 -0.005945 0.007584 0.002214 0.0 0.000469 0.02917 NaN 0.03543 NaN NaN NaN NaN NaN 0.0 -0.0283 0.009003 NaN 0.000377 NaN NaN 0.01635 0.0 NaN 0.0 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN 0.0 NaN NaN NaN NaN 0.006451 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    18.0 -0.0087 0.0 0.0023 0.0153 -0.0055 0.0023 -0.0029 0.004 0.0 NaN 0.0 0.0 0.004 NaN NaN NaN -0.0007 0.0 0.0 0.0 0.0 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN -0.0004 NaN NaN NaN 0.0 0.0 NaN 0.0 NaN -0.00102 NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    19.0 0.0 0.0 0.0 0.0 -0.00878 0.0 0.0 NaN 0.0 NaN 0.0 0.0 -0.002636 0.0097411 0.0 NaN 0.007147 0.0 0.0 NaN 0.0 0.0 NaN 0.0 -0.02978 0.0 NaN 0.0 0.0 NaN 0.0 NaN 0.0 0.0 NaN 0.0 0.0 NaN 0.0 NaN 0.0 0.001118 0.0 NaN 0.0 NaN NaN NaN 0.0 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    20.0 0.010238 0.0 0.0 0.0 0.0 0.0 0.0052371 0.0 0.0 0.0 0.0 NaN 0.0 NaN NaN NaN NaN 0.0 NaN 0.0 0.0 0.0 NaN 0.0 0.0 NaN NaN NaN NaN NaN NaN 0.0 0.0 NaN NaN NaN 0.0 NaN 0.0 NaN 0.0 0.0 0.0 0.0 0.0 NaN 0.0 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    21.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 NaN 0.0 0.0 0.0 NaN 0.001269 NaN NaN NaN -0.003676 0.0 0.0 0.0 0.0 0.0 0.0069046 0.0 0.0 0.0 0.0 0.0 NaN 0.0 NaN NaN 0.0 NaN NaN NaN 0.0 NaN 0.0 NaN NaN 0.0 0.0 0.0 0.0 NaN NaN NaN 0.0 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    22.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 NaN 0.0 0.0 0.0 NaN 0.002745 NaN NaN 0.0 NaN 0.0 0.0 0.0 0.0 0.0 0.0067179 0.0 0.0 0.0 NaN NaN NaN 0.0 NaN 0.0 0.0 NaN 0.0 NaN 0.0 NaN 0.0 NaN NaN 0.0 0.0 NaN 0.0 NaN NaN 0.0 0.0 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    23.0 0.014417 0.0 0.0 0.0 0.0 0.0 NaN NaN 0.0037129 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN 0.0062484 0.003701 0.0 0.0 NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN 0.0 NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    24.0 -0.002416 0.0 -0.003239 0.0 0.0 0.0 0.0 0.0 -0.0025 NaN 0.0 0.0 0.0 0.0 0.0 0.0 0.01466 0.0 0.0 0.0 0.0 0.0 0.0 0.0 -0.03582 0.0 0.0 0.0 NaN 0.0 NaN 0.0 0.0 NaN 0.0 0.0 0.0 0.0 0.0 NaN NaN 0.0 0.0 NaN 0.0 NaN NaN NaN 0.001855 0.0 NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    25.0 0.03317 0.0 0.01661 0.0 0.002987 0.0 0.0 0.0 0.006309 NaN 0.0 0.0 0.01675 0.0 0.0 0.03722 0.0 NaN 0.01158 0.0 0.0 0.0 NaN 0.009219 0.0 0.0 0.0 NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN 0.0 NaN NaN 0.0 0.0 NaN 0.0 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    26.0 -0.006266 0.0 -0.005376 0.004381 0.0 0.001238 0.0 NaN 0.0 NaN 0.0 NaN 0.0 NaN NaN NaN NaN NaN 0.0 NaN 0.0 0.0 NaN 0.0 0.0 0.0 0.0 0.0 NaN NaN NaN 0.0 0.0 0.0 NaN NaN 0.0 NaN 0.0 NaN NaN 0.0 0.0 NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    27.0 0.0 NaN 0.008299 4.51638e-7 0.0 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN 0.0 NaN NaN 0.0 0.0 0.0 0.0 NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    28.0 0.0 0.0 0.003682 -0.01225 NaN 0.01177 NaN NaN 0.0 NaN 0.0 NaN 0.0 NaN NaN NaN NaN NaN 0.0 NaN 0.0 NaN NaN 0.0 NaN 0.0 NaN 0.0 NaN NaN NaN 0.0 NaN NaN NaN NaN 0.0 NaN NaN NaN NaN 0.0 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    29.0 0.0 NaN 0.0 NaN NaN 0.0 NaN NaN 0.0 0.0 NaN 0.0 0.0 0.0 NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN 0.0 NaN NaN NaN 0.0 NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    30.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 NaN 0.0 0.0 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN 0.0 0.0 NaN 0.0 NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    31.0 0.0 NaN -0.00141 0.00717 0.0 0.0 0.0 0.0 0.0 NaN 0.0 NaN NaN NaN NaN NaN 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN 0.0 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    32.0 0.0 NaN 0.0 0.0 0.0 0.0 NaN NaN 0.0 0.0 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN 0.0 NaN 0.0 NaN 0.0 NaN 0.0 NaN 0.0 NaN NaN NaN 0.0 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    33.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 NaN 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 NaN NaN NaN NaN 0.0 0.0 NaN 0.0 NaN NaN NaN NaN NaN 0.0 0.0 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    34.0 0.0 -0.0027 NaN NaN 0.0 NaN NaN NaN 0.0 0.0 NaN NaN 0.0 NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    35.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 NaN 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN 0.0 NaN NaN NaN NaN 0.0 NaN 0.0 NaN 0.0 NaN 0.0 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN 0.0 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    36.0 0.0 0.0 NaN NaN 0.0 0.0 0.0 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    37.0 0.0 0.0 0.007813 0.0 0.0 0.0 NaN NaN 0.0 0.0 0.0 0.0 0.0 NaN NaN NaN NaN NaN 0.0 0.0 0.0 0.0 0.0 0.0 NaN 0.0 NaN 0.0 NaN 0.0 NaN NaN NaN NaN NaN 0.0 0.0 NaN 0.0 NaN 0.0 0.0 0.013871 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    38.0 0.0 0.0 -0.007754 0.0 0.0 0.0 NaN NaN NaN NaN NaN NaN 0.0 NaN 0.0 0.0 NaN 0.0027 NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 0.0 0.0 NaN 0.0 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    39.0 0.0 0.0 0.000668 0.0 0.0 0.0 0.0 NaN 0.011386 0.0 0.0 0.0 NaN 0.0 0.0 0.0 0.0 NaN 0.0 0.0 0.0 0.0 NaN 0.0 0.0 0.0 NaN NaN 0.0 NaN 0.0 NaN NaN 0.0 0.0 0.0 0.0 0.0 0.0 NaN NaN 0.0 0.0 0.0 NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    40.0 0.0 NaN 0.0 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    41.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 NaN 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN 0.0 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN 0.0 NaN NaN NaN 0.0 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    42.0 0.001291 0.0 -0.001557 0.0 -0.006309 0.0 0.003455 0.0 -0.004653 0.0 4.1e-5 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.004586 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 NaN NaN 0.0 NaN NaN NaN 0.0 0.0 0.0 0.0 0.0 0.0 0.0 NaN 0.0 0.0 0.0 0.0 0.0 0.0 NaN 0.0 0.0 0.0 0.0 0.0119213 0.0 0.0076136 0.0 0.0 0.0 0.0360498 NaN
    43.0 0.0024016 0.0037455 0.0 0.0 -0.008186501 0.0053598 -7.437899e-14 NaN 0.0 NaN 0.0 NaN 0.0 0.0 0.0 NaN NaN 0.0 0.0 0.0 0.0 0.0 NaN 0.0 0.0 0.0 NaN 0.0 NaN NaN NaN NaN 0.0 NaN NaN NaN 0.00091722 0.0 0.0 NaN NaN 0.0 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN 0.0 NaN 0.0 0.0 NaN -0.0080098 0.0 0.0 NaN NaN NaN
    44.0 0.0 NaN -0.015439 0.0 NaN NaN 0.0 NaN NaN NaN 0.0 0.0 0.0 NaN NaN NaN NaN NaN NaN 0.0 0.0 NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    45.0 0.0 0.0 0.0 0.0 0.0 0.025981 0.0 NaN 0.0045351 NaN 0.0 0.0 0.0 NaN NaN 0.0 NaN 0.0 0.0 0.0 0.0 0.0 NaN 0.0 0.0 NaN NaN 0.0 NaN NaN NaN 0.0 NaN NaN 0.0 NaN 0.0 NaN NaN NaN NaN 0.0 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    46.0 0.0 0.0 0.0 0.0 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN 0.0 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    47.0 0.0 0.0 0.001331 0.0 0.004691 0.0 -0.001586 NaN 0.0 NaN 0.0 NaN NaN 0.0 0.0 0.0 -0.009617 -0.00015 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN 0.0 0.0 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    48.0 -0.000623 0.0 2.0e-6 -0.000137 -0.005877 -0.00208 0.0 NaN 0.0 NaN 0.0 NaN 0.0 NaN NaN 0.0 NaN NaN NaN 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN 0.0 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    49.0 0.000583 0.0 0.000281 0.000607 0.0 0.0 -0.000916 NaN 0.0 NaN 0.0 NaN -4.0e-6 NaN NaN NaN NaN NaN 0.0 NaN 0.0 0.0 NaN -0.001271 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN 0.0 NaN NaN 0.0 NaN NaN NaN NaN 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    52.0 0.0 0.0 0.0 0.0 -0.035 0.0 NaN 0.0 0.0 NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN 0.0 NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    53.0 0.0 0.0 0.0 0.0 -0.00244 -0.00166 -0.0048 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    55.0 4.0e-5 0.0 -0.00782 0.00633 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 -6.5e-5 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN
    56.0 0.0 0.0 0.00242 0.01424 0.001188 0.01159 0.0 NaN NaN NaN 0.0126 NaN 0.00839 0.0 0.008006 -0.00165 NaN NaN 0.0 NaN 0.0 0.0 0.0 0.0 NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN 0.0 0.0 NaN 0.0 NaN NaN NaN NaN NaN 0.0 -0.000969 0.0 -0.0016 NaN NaN NaN NaN NaN NaN NaN NaN NaN
    61.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 NaN NaN NaN 0.013472 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0078 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN
    84.0 0.0 0.0 0.0 -0.0007792 0.0 -0.0098745 0.0019814 NaN 0.0047161 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 0.0094424 NaN 0.0005537 NaN -0.0025991 NaN 0.0 -0.0518203
    85.0 0.0095873 0.0 0.0025793 0.0 0.0 0.0 0.0 NaN 0.0074865 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0018563 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0208006 0.0 0.0382461 NaN 0.0 NaN NaN NaN NaN
    87.0 -0.0055872 0.0 -0.0027586 -0.0186009 0.0 0.0 NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0160959 0.0 NaN NaN 0.0 NaN NaN NaN
    89.0 0.0010454 -0.0254894 0.0 0.0 -0.0134709 -0.0062036 -0.0067585 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0137839 -0.0351738 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 2.63e-5 NaN NaN 0.0 -0.0092342 NaN NaN NaN NaN
    90.0 -0.0081757 0.0 0.0 0.0 0.0 0.0 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN -0.0356672 0.0 NaN NaN NaN NaN
    91.0 0.0 0.0056279 -0.0001404 0.0 0.0 0.0 0.0 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -0.0035024 NaN 0.0 NaN NaN 0.0 NaN NaN NaN
    93.0 0.0 0.0 0.0 0.0 0.0 0.0 0.0 NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN
    98.0 0.0 0.0 0.0 0.0 0.000518 0.0 -0.000807 NaN -0.0082641 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -0.029654 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN 0.0 NaN
    99.0 -0.0330202 -0.0078001 0.0028389 -0.0217293 -0.0148923 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -0.0073933 NaN NaN NaN NaN NaN NaN NaN 0.0
]

Rk_arr_glob = [
    0.6325
    0.6325
    0.6325
    0.6325
    1.2832
    1.2832
    1.2832
    1.2832
    1.2832
    0.3763
    0.3763
    0.91
    0.91
    0.91
    1.2302
    1.063
    0.6895
    0.8585
    1.7334
    1.08
    1.7048
    1.7048
    0.7173
    1.27
    1.27
    1.9
    1.1434
    1.1434
    1.1434
    1.6607
    1.6607
    1.6607
    1.6607
    1.368
    1.368
    1.368
    1.0746
    1.0746
    1.1849
    1.4578
    1.2393
    1.0731
    1.5575
    1.5575
    0.8
    0.9919
    0.9919
    0.9919
    1.8
    1.8
    1.8
    2.65
    2.618
    0.5365
    2.644
    2.5
    2.887
    0.4656
    1.24
    1.289
    1.535
    1.299
    2.088
    1.076
    1.209
    0.9214
    1.303
    3.6
    1
    0.5229
    0.8814
    2
    2.381
    1.284
    1.284
    0.8215
    1.6
    0.7136
    0.3479
    0.347
    1.7023
    1.4046
    1.0413
    0.8
    2.45
    3.981
    3.7543
    3.5268
    3.2994
    1.4515
    1.5
    1.5
    2.4748
    2.2739
    2.0767
    2.4617
    2.4617
    1.7943
    1.6282
    1.4621
    1.3601
    0.683
    0.9104
    1.063
    0.9104
    2.42
    2.42
    2.42
    2.687
    2.46
    1.613
    1.3863
    1.1589
    1.3662
    1.843
    5.621
    2.7867
    3.9628
    2.1094
    2.4873
    3.371
    1.0678
    0.9903
    1.5654
    3.8183
]

Qk_arr_glob = [
    1.0608
    0.7081
    0.3554
    0
    1.6016
    1.2489
    1.2489
    0.8962
    0.4582
    0.4321
    0.2113
    0.949
    0.7962
    0.3769
    0.8927
    0.8663
    0.8345
    0.9938
    2.4561
    0.975
    1.67
    1.5542
    0.771
    1.6286
    1.4228
    1.8
    1.6022
    1.2495
    0.8968
    1.6904
    1.3377
    0.985
    0.985
    1.4332
    1.0805
    0.7278
    1.176
    0.824
    0.8067
    0.9022
    0.633
    0.353
    1.5193
    1.1666
    0.9215
    1.3654
    1.0127
    0.66
    2.5
    2.1473
    1.7946
    2.3778
    3.1836
    0.3177
    2.5
    2.304
    2.241
    0.3589
    1.068
    1.762
    1.316
    1.289
    2.4
    0.9169
    1.4
    1.3
    1.132
    2.692
    0.92
    0.7391
    0.7269
    2.093
    1.522
    1.266
    1.098
    0.5135
    0.9
    0.8635
    0.1071
    0
    1.8784
    1.4
    1.0116
    1.2742
    2.8912
    3.2
    2.892
    2.58
    2.352
    1.248
    1.08
    1.08
    1.9643
    1.5754
    1.1866
    2.192
    1.842
    1.34
    1.06
    0.78
    1.8031
    0.3418
    0.6538
    1.123
    0.6538
    2.4976
    2.0018
    2.2497
    2.12
    1.808
    1.368
    1.06
    0.748
    0.6797
    1.6997
    5.9463
    2.7723
    0.6214
    2.5106
    2.4457
    2.0001
    2.244
    3.5249
    3.8076
    3.6018
]

Idxarr_glob = [
    1.0,
    1.0,
    1.0,
    1.0,
    2.0,
    2.0,
    2.0,
    2.0,
    2.0,
    3.0,
    3.0,
    4.0,
    4.0,
    4.0,
    5.0,
    5.0,
    5.0,
    6.0,
    7.0,
    8.0,
    9.0,
    9.0,
    10.0,
    11.0,
    11.0,
    12.0,
    13.0,
    13.0,
    13.0,
    14.0,
    14.0,
    14.0,
    14.0,
    15.0,
    15.0,
    15.0,
    16.0,
    16.0,
    17.0,
    18.0,
    18.0,
    18.0,
    19.0,
    19.0,
    20.0,
    21.0,
    21.0,
    21.0,
    22.0,
    22.0,
    22.0,
    23.0,
    24.0,
    25.0,
    26.0,
    26.0,
    26.0,
    27.0,
    28.0,
    29.0,
    29.0,
    30.0,
    31.0,
    32.0,
    33.0,
    34.0,
    34.0,
    35.0,
    36.0,
    37.0,
    38.0,
    39.0,
    39.0,
    40.0,
    40.0,
    40.0,
    41.0,
    42.0,
    42.0,
    42.0,
    43.0,
    43.0,
    43.0,
    44.0,
    45.0,
    46.0,
    46.0,
    46.0,
    46.0,
    47.0,
    47.0,
    47.0,
    48.0,
    48.0,
    48.0,
    49.0,
    49.0,
    50.0,
    50.0,
    50.0,
    51.0,
    51.0,
    51.0,
    51.0,
    51.0,
    52.0,
    52.0,
    52.0,
    53.0,
    53.0,
    54.0,
    54.0,
    54.0,
    55.0,
    55.0,
    56.0,
    57.0,
    58.0,
    59.0,
    59.0,
    60.0,
    61.0,
    62.0,
    62.0,
    63.0,
]


end
