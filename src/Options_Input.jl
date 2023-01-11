
#include("Fugacity_modiles/Fugacity_func_list.jl")

mutable struct opts
    
    activity_mdl::String
    fugacity_mdl::String
    mixing::String
    poynting::Bool

end

mutable struct component
    Name::String
    UNIFAC_comp::Matrix{Int64}
    UNIFACmod_comp::Matrix{Int64}
    Antoine_coeff::Vector{Float64}
    Tc::Float64
    Pc::Float64
    ω::Float64
    v_l::Float64
end


function def_component(;Name="Water",UNIFAC_comp=[17 1], UNIFACmod_comp=[19 1], Antoine_coeff=[4.65,1435.26,-64.848],Tc=647.0,Pc=2.212e+7,ω=0.334,v_l=1.8798e-5)

    return component(Name,UNIFAC_comp,UNIFACmod_comp,Antoine_coeff,Tc,Pc,ω,v_l)
    
end


function def_calc_options(;activity_mdl="unity",fugacity_mdl="unity",mixing_mdl="none",poynting=false)
    
    activity_mdl=lowercase(activity_mdl)
    fugacity_mdl=lowercase(fugacity_mdl)
    mixing_mdl=lowercase(mixing_mdl)

    activity_mdl_opts=["unity","unifac","unifacmod"]
    fugacity_mdl_opts=["unity","vdw","rk","srk","pr"]
    mixing_mdl_opts=["none","vdw","k-matrix"]
    poynting_opts=[true, false]
    
    if any(activity_mdl.==activity_mdl_opts)==false
        println(string(activity_mdl," is not a valid option for activity_mdl, reverting to default."))
        activity_mdl="unity"
    elseif any(fugacity_mdl.==fugacity_mdl_opts)==false
        println(string(fugacity_mdl," is not a valid option for fugacity_mdl, reverting to default."))
        fugacity_mdl="unity"
    elseif any(mixing_mdl.==mixing_mdl_opts)==false
        println(string(mixing_mdl," is not a valid option for mixing_mdl, reverting to default."))
        mixing_mdl="none"
    elseif any(poynting.==poynting_opts)==false
        println(string(poynting," is not a valid option for poynting, reverting to default."))
        poynting=false
    end

    return opts(activity_mdl,fugacity_mdl,mixing_mdl,poynting)
end


#UNIFAC function with struct as input
function UNIFAC(struct_vec::Vector{component}, T_k::Float64, x_arr::Vector{Float64})

    M_local = Vector{Matrix{Int64}}(undef, length(struct_vec))
    n_data = 0
    for i in eachindex(struct_vec)
        M_local[i] = struct_vec[i].UNIFAC_comp
        if length(M_local[i]) == 0
            n_data = 1
            println("No data in field UNIFAC_comp for component", struct_vec[i].Name, ". Cannot perform calculation")
        end
    end

    if n_data == 0
        return Activity_UNIFAC.UNIFAC(T_k, M_local, x_arr)
    else
        println("One or more pieces of data necessary to calculation are missing")
    end
end


function UNIFACmod(struct_vec::Vector{component}, T_k::Float64, x_arr::Vector{Float64})

    M_local = Vector{Matrix{Int64}}(undef, length(struct_vec))
    n_data = 0
    for i in eachindex(struct_vec)
        M_local[i] = struct_vec[i].UNIFACmod_comp
        if length(M_local[i]) == 0
            n_data = 1
            println("No data in field UNIFACmod_comp for component", struct_vec[i].Name, ". Cannot perform calculation")
        end
    end

    if n_data == 0
        return Activity_UNIFACmod.UNIFACmod(T_k, M_local, x_arr)
    else
        println("One or more pieces of data necessary to calculation are missing")
    end
end

function EOS_Pure(struct_in::component;T=298.15,P=101.3e+3,model="SRK")
   

    if model=="Ideal"
        res=(1,1)
    elseif model=="VdW"
        res=Fugacity.VdW_EOS(P,T,struct_in.Pc,struct_in.Tc)
    elseif model=="RK"
        res=Fugacity.RK_EOS(P,T,struct_in.Pc,struct_in.Tc)
    elseif model=="SRK"
        res=Fugacity.SRK_EOS(P,T,struct_in.Pc,struct_in.Tc,struct_in.ω)
    elseif model=="PR"
        res=Fugacity.PR_EOS(P,T,struct_in.Pc,struct_in.Tc,struct_in.ω)
    end

    return res
end
    


