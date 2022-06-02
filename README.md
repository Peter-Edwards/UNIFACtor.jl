# UNIFACTOR.jl

UNIFACTOR.jl is a package for Julia that allows a used to easily calculate the activity coefficients of an arbitrary number of different components of a liquid mixture.

The current iteration of the package can only use the original UNIFAC method to find the activity coefficients. Any and all of the UNIFAC coefficients used in the model were taken from the official list published by the UNIFAC Consortium: 

http://unifac.ddbst.de/published-parameters-unifac.html

the first 50 groups UNIFAC groups are included in the package’s UNIFAC model along with 108 different functional groups.

## Usage

Currently the only available function is called using:

```julia
UNIFAC.Activity(T_k,M_lst,x_arr)
```

T_k is a temperature value in kelvin. The x_arr is a 1xn matrix of molar faction values (where n is the number of components). M_lst is a vector of 1xm matrices describing the composition of each of the component molecules, where m is dependant on the functional groups present. As stated previously, 108 different functional groups are available for use, however the matrices present in in M_lst only need to be long enough to fully describe the molecule. An example of this is given below:

```julia
x1=0.5
x2=0.25
x3=1-x1-x2
x_arr=[x1 x2 x3]

T_k=298

M_1=[0 0 0 0 0 0 0 0 0 6]; #benzene
M_2=[1 1 0 0 0 0 0 0 0 0 0 0 0 0 1]; #ethanol
M_3=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1]; #acetone

M_lst=[M_1,M_2,M_3]

UNIFAC.Activity(T_k,M_lst,x_arr)
```
In each case it is assumed that all of the remaining functions groups out of the 108 are treated as zero.

It should also be noted that the indexing of functional groups differs slightly from the ddbst website in that the C=C groups is not indexed as 70, but instead indexed as 9. the order of molecules is as follows:

CH3,CH2,CH,C,CH2=CH,CH=CH,CH2=C,CH=C,C=C,ACH,AC,ACCH3,ACCH2,ACCH,OH(alcohol),CH3OH (Methanol),H2O,ACOH (Alchohol),CH3CO (ketone),CH2CO(Ketone), CHO (Aldheyde),...

it may be observed that after C=C the order is exactly the same as listed on the ddbst website.
