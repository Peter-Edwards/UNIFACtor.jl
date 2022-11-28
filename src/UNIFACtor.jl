module UNIFACtor

#Activity functions
include("Activity_modules/Activity_UNIFAC.jl")
include("Activity_modules/Activity_UNIFACmod.jl")
using .Activity_UNIFACmod, .Activity_UNIFAC
export UNIFACmod, UNIFAC
#Fugacity Functions
#include("Fugacity_modiles/Fugacity_func_list.jl")

#system parameter functions


end

