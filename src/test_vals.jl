
using UNIFACtor

M_U = [[17 1], [1 1; 2 1; 15 1]]
M_UM = [[19 1], [1 1; 2 1; 15 1]]
x_v = [0.128, 0.872]
y_v = x_v
T_k = 351.35
P = 101.3e+3
Pc = [221.2 * 100e+3, 63 * 100e+3]
Tc = [647, 513.9]
ω = [0.344, 0.644]
A_coeff = [[4.65, 1435.26, -64.848], [5.24677, 1598.673, -46.424]]
A_coeff3 = [[4.65, 1435.26, -64.848], [5.24677, 1598.673, -46.424], [4.42448, 1312.253, -32.445]]

#setting up component structs
water_struct = def_component()
eth_struct = def_component(Name="Ethanol", UNIFAC_comp=M_U[2], UNIFACmod_comp=M_UM[2], Pc=Pc[2], Tc=Tc[2], Antoine_coeff=A_coeff[2], ω=ω[2])
ace_struct = def_component(Name="Acetone", UNIFAC_comp=[1 1; 19 1], UNIFACmod_comp=[1 1; 21 1], Pc=46.9 * 100e+3, Tc=508.1, ω=0.304, Antoine_coeff=[4.42448, 1312.253, -32.445])
benz_struct = def_component(Name="Benzene", UNIFAC_comp=[10 6], UNIFACmod_comp=[10 6], Pc=4.89e+6, Tc=562.1, ω=0.212, Antoine_coeff=[4.72583, 1660.652, -1.461])
hex_struct = def_component(Name="Hexane", UNIFAC_comp=[1 2; 2 4], UNIFACmod_comp=[1 2; 2 4], Pc=3.01e+6, Tc=507, ω=0.299, Antoine_coeff=[4.00266, 1171.53, -48.784])

comp_struct = [water_struct, eth_struct]
#comp_struct3=[water_struct,eth_struct,ace_struct]
comp_struct3 = [ace_struct, benz_struct, water_struct]

calc_options = def_calc_options(; activity_mdl="unifacmod", fugacity_mdl="srk")
x_arr=[0.715,0.181]
comp_tst=[ace_struct,hex_struct,water_struct]
