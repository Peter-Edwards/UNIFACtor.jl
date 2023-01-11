

module Activity_UNIFAC

export UNIFAC
using DelimitedFiles

#the overall activity function
function UNIFAC(T_k::Float64, M_lst::Vector{Matrix{Int64}}, x_arr::Vector{Float64})

    #normalises x_arr
    x_arr = x_arr ./ sum(x_arr)
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

    return (phi, theta)
end

#finding the overall combinatorial contributions
function func_comb_mk2(x_arr, M_local)
    #(Ri, Qi) = func_qiri_mk2(M_local,Qk_arr_glob,Rk_arr_glob)
    (Ri, Qi) = func_qiri_mk3(M_local)
    (phi, theta) = func_thetaphi(Ri, Qi, x_arr)

    ln_gamc = Vector{Float64}(undef, length(Ri))

    @. ln_gamc = 1 - phi + log(phi) - 5 * Qi * (1 - phi / theta + log(phi / theta))
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

    # for i in axes(M_local_comp,1)
    #     sum_frac[i]/sum_frac_comp
    # end
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
                Anm[info_in.Group[i], info_in.Group[j]]
                psi_mat[i, j] = exp(-(Anm[info_in.Group[i], info_in.Group[j]] / T))
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

#Group Interraction Parameters (50x50)
data_dir=joinpath(@__DIR__, "Parameter_data")
#Anm = readdlm("Parameter_data/UNIFAC_Int_Params.txt")
Anm = readdlm(joinpath(data_dir,"UNIFAC_Int_Params.txt"))
#structure params (Rk,Qk,group index)
#structure_params = readdlm("Parameter_data/UNIFAC_Grp_Params.txt")
structure_params = readdlm(joinpath(data_dir,"UNIFAC_Grp_Params.txt"))
Qk_arr_glob = structure_params[1, :]
Rk_arr_glob = structure_params[2, :]
Idxarr_glob = structure_params[3, :]


end
