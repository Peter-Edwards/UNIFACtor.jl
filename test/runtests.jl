using Test
using UNIFACtor

x1 = 0.5
x2 = 0.25
x3 = 1 - x1 - x2
x_arr = [x1, x2, x3]

T_k = 298.5

M_1 = [10 6] #benzene
M_2 = [1 1; 2 1; 15 1]#ethanol
M_3 = [1 1; 19 1]#actone
M_3D = [1 1; 21 2]#dortmund acetone
#M_W=[17 1]
#M_WD=[19 1] #dortmund water

M_lst = [M_1, M_2, M_3]
test_UNIFAC1 = UNIFAC(T_k, M_lst, x_arr)

M_lst = [M_1, M_2, M_3D]
test_UNIFACmod1 = UNIFACmod(T_k, M_lst, x_arr)

testmean1 = sum((test_UNIFAC1 - test_UNIFACmod1) .^ 2 ./ test_UNIFAC1) / length(test_UNIFAC1)

println(test_UNIFAC1, test_UNIFACmod1)
@test testmean1 < 0.05

