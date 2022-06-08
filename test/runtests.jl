using Test
using UNIFACtor

x1=0.5
x2=0.25
x3=1-x1-x2
x_arr=[x1 x2 x3]

T_k=298.5

M_1=[0 0 0 0 0 0 0 0 0 6]; #benzene
M_2=[1 1 0 0 0 0 0 0 0 0 0 0 0 0 1]; #ethanol
M_3=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]; #acetone
M_3D=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]; #acetone dortmund

M_lst=[M_1,M_2,M_3]
test_UNIFAC1=UNIFAC.Activity(T_k,M_lst,x_arr)

M_lst=[M_1,M_2,M_3D]
test_UNIFACmod1=UNIFACmod.Activity(T_k,M_lst,x_arr)

testmean1=sum((test_UNIFAC1-test_UNIFACmod1).^2 ./test_UNIFAC1)/length(test_UNIFAC1)


@test testmean1<0.05