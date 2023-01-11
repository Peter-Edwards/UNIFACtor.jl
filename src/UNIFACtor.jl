module UNIFACtor

#Activity functions
include("Activity_modules/Activity_UNIFAC.jl")
include("Activity_modules/Activity_UNIFACmod.jl")
#Fugacity Functions
include("Fugacity_modiles/Fugacity_func_list.jl")

#Input and Options functions
export UNIFAC
export UNIFACmod
export def_component
export def_calc_options
export EOS_Pure
include("Options_Input.jl")
include("EqmCalc_func_list.jl")
#system parameter functions


end

