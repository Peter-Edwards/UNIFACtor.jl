

module Input

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
end


function def_component(;Name="Water",UNIFAC_comp=[17 1], UNIFACmod_comp=[19 1], Antoine_coeff=[4.65,1435.26,-64.848],Tc=647.0,Pc=2.212e+7,ω=0.334)

    return component(Name,UNIFAC_comp,UNIFACmod_comp,Antoine_coeff,Tc,Pc,ω)
    
end


function def_calc_options(;activity_mdl="unity",fugacity_mdl="unity",mixing_mdl="none",poynting=false)

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




end
