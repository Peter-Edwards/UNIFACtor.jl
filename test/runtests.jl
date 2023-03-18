using Test
using UNIFACtor

x_arr = [0.5, 0.25, 0.25]
T_k = 298.5

#M_1 = [10 6] #benzene
#M_2 = [1 1; 2 1; 15 1]#ethanol
#M_3 = [1 1; 19 1]#actone
#M_3D = [1 1; 21 1]#dortmund acetone
#M_W=[17 1]
#M_WD=[19 1] #dortmund water
M_1 = def_component(UNIFAC_comp=[10 6], UNIFACmod_comp=[10 6])
M_2 = def_component(UNIFAC_comp=[1 1; 2 1; 15 1], UNIFACmod_comp=[1 1; 2 1; 15 1])
M_3 = def_component(UNIFAC_comp=[1 1; 19 1], UNIFACmod_comp=[1 1; 21 1])

M_lst = [M_1, M_2, M_3]
test_UNIFAC1 = UNIFAC(M_lst, T_k, x_arr)
test_UNIFACmod1 = UNIFACmod(M_lst, T_k, x_arr)

testmean1 = sum((test_UNIFAC1 - test_UNIFACmod1) .^ 2 ./ test_UNIFAC1) / length(test_UNIFAC1)

println(test_UNIFAC1, test_UNIFACmod1)
@test testmean1 < 0.005
@test isapprox(test_UNIFAC1, [1.32, 2.15, 1.04]; atol=0.1)
@test isapprox(test_UNIFACmod1, [1.32, 2.15, 1.04]; atol=0.1)
