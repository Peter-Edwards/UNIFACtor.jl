
module UNIFAC
#the overall activity function

function Activity(T_k,M_lst,x_arr)
#list includes groups as follows (in order)
#CH3,CH2,CH,C,CH2=CH,CH=CH,CH2=C,CH=C,C=C,ACH,AC,ACCH3,ACCH2,ACCH,OH(alcohol),CH3OH (Methanol),H2O,ACOH (Alchohol),CH3CO (ketone),CH2CO(Ketone), CHO (Aldheyde),...
M_lst=func_M_lst(M_lst)


#changes zero values so the the math is doable. useful for generating PXY graphs and such.
x_arr[x_arr .<= 1e-8] .= 1e-8

#residual contribution
#forming a matrix of group information for mixture and pure substance interraction
grp_tab=func_infotab_grpV2(f_Qk(),M_lst,x_arr)
ind_tab=func_infotab_indV2(f_Qk(),M_lst)


#updating each struct with the LnRho line 
grp_tab=func_lnRho(grp_tab,T_k)
for i in eachindex(ind_tab)
ind_tab[i]=func_lnRho(ind_tab[i],T_k)
end

idx=Array{Float64}(undef,0)
for i in eachindex(ind_tab)
    push!(idx,func_lngamV2(ind_tab[i],grp_tab))
end

#combining combinatorial and residual contributions
gam_arr=exp.(vec(func_comb_mk2(x_arr,M_lst))+idx)

return gam_arr

end



#finding Ri and Qi for every molecule
function func_qiri_mk2(M_lst::Matrix,Qk::Array,Rk::Array)

    Ri=sum(Rk.*M_lst', dims=1)
    Qi=sum(Qk.*M_lst', dims=1)
    return (Ri,Qi)

end

#finding theta and phi
function func_thetaphi(Ri::Array,Qi::Array,x_arr::Array)
    phi=Ri./sum(Ri.*x_arr)
    theta=Qi./sum(Qi.*x_arr)

    return (phi,theta)
end


function func_comb_mk2(x_arr,M_lst)
    (Ri,Qi)=func_qiri_mk2(M_lst,f_Qk(),f_Rk())
    (phi,theta)=func_thetaphi(Ri,Qi,x_arr)

    ln_gamc=zeros(1,length(Qi))

    @. ln_gamc=1-phi+log(phi)-5*Qi*(1-phi/theta+log(phi/theta))
    return ln_gamc
end

#the func_idxarr function is just used to generate the group indexes for f_idxarr

# function func_idxarr()
#     size_arr=[4 5 2 3 1 1 1 1 2 1 2 1 4 3 3 2 1 3 2 2 3 3 2 1 1 3 1 1 2 1 1 1 1 2 1 1 1 1 2 3 1 4 3 1 8 6 2 3 1 3]
#     idx_arr=Array{Float64}(undef,0)
#     ind=1
#     for i in 1:length(size_arr)
#         append!(idx_arr,ones(size_arr[i])*ind) 
#         ind=ind+1
#     end
#     return idx_arr
# end

#function to assign group to each molecule
function func_assign_grp(M_lst,idx_arr)
    (n,m)=size(M_lst)

    grp_assign=zeros(n,m)
    for i in 1:n
        M_fcus=M_lst[i,:]
        grp_assign[i,:]=map(!=(0),M_fcus).*idx_arr
    end
    return grp_assign
end 

#find area fraction for each pure molecule and for the mixture
function func_arrfrac(Qk,M_lst,x_arr)
    M_lst=M_lst.*x_arr'
    M_ex=sum(M_lst, dims=1)
    M_ex=M_ex/sum(M_ex)
    Ar_frac=M_ex.*Qk'
    Ar_denom=sum(Ar_frac)
    Ar_frac=Ar_frac./Ar_denom

    return Ar_frac


end

#the struct for infotab
struct InfoStruct

    index::Array{Int64}
    Quan::Array{Int64}
    Group::Array{Int64}
    AreaFrac::Array{Float64}
    Qk::Array{Float64}
    LnRho::Array{Float64}

end
#building a matrix where the information for each pure component is held. I want to change the code so this matrix turns into a struct, it just seems more elegant and less of a bodge. 
#thats what 4 years of mainly matlab will do to a man.

#new and improved version of infotab_ind (with structs)
function func_infotab_indV2(Qk,M_lst)
    idx_s=f_idxarr()
    tab_out=Array{InfoStruct}(undef,0)
    for i in axes(M_lst)[1]
        fcus_Mlst=M_lst[i,:]'
        index=findall(x->x>=1,vec(fcus_Mlst))
        Arr_frac=func_arrfrac(Qk,fcus_Mlst,1)

        # list in decending order, index of group, quantity of group, interaction catagory of group, Area Fraction of group, Qk value of group
        tab_fcus=InfoStruct(index,fcus_Mlst[index],idx_s[index],Arr_frac[index],Qk[index],[])
        append!(tab_out,[tab_fcus])
    end
    return tab_out

end

#building a struct where the information for all of the group mixture is held
function func_infotab_grpV2(Qk,M_lst,x_arr)

    Arr_frac=func_arrfrac(Qk,M_lst,x_arr)
    idx_gs=f_idxarr()

    Quan_idx=sum((M_lst),dims=1)
    index=findall(x->x>=1,vec(Quan_idx))
    tab_out=InfoStruct(index,Quan_idx[index],idx_gs[index],Arr_frac[index],Qk[index],[])
    # list in decending order, index of group, quantity of group, interaction catagory of group, Area Fraction of group, Qk value of group
    
    return tab_out
end



#a function to find psi using the interraction parameter matrix and system temperature
function func_intparams(T,n,m)

Int_GP=[0.0 86.02 61.13 76.5 986.5 697.2 1318.0 1333.0 476.4 677.0 232.1 507.0 251.5 391.5 255.7 206.6 920.7 287.77 597.0 663.5 35.93 53.76 24.9 104.3 11.44 661.5 543.0 153.6 184.4 354.55 3025.0 335.8 479.5 298.9 526.5 689.0 -4.189 125.8 485.3 -2.859 387.1 -450.4 252.7 220.3 -5.869 390.9 553.3 187.0 216.1 92.99; -35.36 0.0 38.81 74.15 524.1 787.6 270.6 526.1 182.6 448.8 37.85 333.5 214.5 240.9 163.9 61.11 749.3 280.5 336.9 318.9 -36.87 58.55 -13.99 -109.7 100.1 357.5 NaN 76.302 NaN 262.9 NaN NaN 183.8 31.14 179.0 -52.87 -66.46 359.3 -70.45 449.4 48.33 NaN NaN 86.46 NaN 200.2 268.1 -617.0 62.56 NaN; -11.12 3.446 0.0 167.0 636.1 637.35 903.8 1329.0 25.77 347.3 5.994 287.1 32.14 161.7 122.8 90.49 648.2 -4.449 212.5 537.4 -18.81 -144.4 -231.9 3.0 187.0 168.0 194.9 52.07 -10.43 -64.69 210.4 113.3 261.3 154.26 169.9 383.9 -259.1 389.3 245.6 22.67 103.5 -432.3 238.9 30.04 -88.11 NaN 333.3 NaN -59.58 -39.16; -69.7 -113.6 -146.8 0.0 803.2 603.25 5695.0 884.9 -52.1 586.6 5688.0 197.8 213.1 19.02 -49.29 23.5 664.2 52.8 6096.0 872.3 -114.1 -111.0 -80.25 -141.3 -211.0 3629.0 4448.0 -9.451 393.6 48.49 4975.0 259.0 210.0 -152.55 4284.0 -119.2 -282.5 101.4 5629.0 -245.39 69.26 683.3 355.5 46.38 NaN NaN 421.9 NaN -203.6 184.9; 156.4 457.0 89.6 25.82 0.0 -137.1 353.5 -259.7 84.0 -203.6 101.1 267.8 28.06 83.02 42.7 -323.0 -52.39 170.0 6.712 199.0 75.62 65.28 -98.12 143.1 123.5 256.5 157.1 488.9 147.5 -120.5 -318.9 313.5 202.1 727.8 -202.1 74.27 225.8 44.78 -143.9 NaN 190.3 -817.7 202.7 -504.2 72.96 -382.7 -248.3 NaN 104.7 57.65; 16.51 -12.52 -50.0 -44.5 249.1 0.0 -181.0 -101.7 23.39 306.4 -10.72 179.7 -128.6 359.3 -20.98 53.9 489.7 580.5 53.28 -202.0 -38.32 -102.5 -139.4 -44.76 -28.25 75.14 457.88 -31.09 17.5 -61.76 -119.2 212.1 106.3 -119.1 -399.3 -5.224 33.47 -48.25 -172.4 NaN 165.7 NaN NaN NaN -52.1 NaN NaN 37.63 -59.4 -46.01; 300.0 496.1 362.3 377.6 -229.1 289.6 0.0 324.5 -195.4 -116.0 72.87 233.87 540.5 48.89 168.0 304.0 459.0 459.0 112.6 -14.09 325.4 370.4 353.7 497.5 133.9 220.6 399.5 887.1 NaN 188.0 12.72 NaN 777.1 NaN -139.0 160.8 NaN NaN 319.0 NaN -197.5 -363.8 NaN -452.2 NaN 835.6 139.6 NaN 407.9 NaN; 275.8 217.5 25.34 244.2 -451.6 -265.2 -601.8 0.0 -356.1 -271.1 -449.4 -32.52 -162.9 -832.97 NaN NaN -305.5 -305.5 NaN 408.9 NaN 517.27 NaN 1827.0 6915.0 NaN -413.48 8484.0 NaN NaN -687.1 NaN NaN NaN NaN NaN NaN NaN NaN NaN -494.2 NaN NaN -659.0 NaN NaN NaN NaN NaN 1005.0; 26.76 42.92 140.1 365.8 164.5 108.7 472.5 -133.1 0.0 -37.36 -213.7 -190.4 -103.6 NaN -174.2 -169.0 6201.0 7.341 481.7 669.4 -191.7 -130.3 -354.6 -39.2 -119.8 137.5 548.5 216.1 -46.28 -163.7 71.46 53.59 245.2 -246.6 -44.58 -63.5 -34.57 NaN -61.7 NaN -18.8 -588.9 NaN NaN NaN NaN 37.54 NaN NaN -162.6; 505.7 56.3 23.39 106.0 529.0 -340.2 480.8 -155.6 128.0 0.0 -110.3 766.0 304.1 NaN NaN NaN NaN NaN -106.4 497.5 751.9 67.52 -483.7 NaN NaN NaN NaN NaN NaN NaN NaN 117.0 NaN 2.21 NaN -339.2 172.4 NaN -268.8 NaN -275.5 NaN NaN NaN NaN NaN NaN NaN NaN NaN; 114.8 132.1 85.84 -170.0 245.4 249.63 200.8 -36.72 372.2 185.1 0.0 -241.8 -235.7 NaN -73.5 -196.7 475.5 -0.13 494.6 660.2 -34.74 108.9 -209.7 54.57 442.4 -81.13 NaN 183.0 NaN 202.3 -101.7 148.3 18.88 71.48 52.08 -28.61 -275.2 NaN 85.33 NaN 560.2 NaN NaN NaN NaN NaN 151.8 NaN NaN NaN; 329.3 110.4 18.12 428.0 139.4 227.8 124.63 -234.25 385.4 -236.5 1167.0 0.0 -234.0 NaN NaN NaN NaN -233.4 -47.25 -268.1 NaN 31.0 -126.2 179.7 24.28 NaN NaN NaN 103.9 NaN NaN NaN 298.13 NaN NaN NaN -11.4 NaN 308.9 NaN -70.24 NaN NaN NaN NaN NaN NaN NaN NaN NaN; 83.36 26.51 52.13 65.69 237.7 238.4 -314.7 -178.5 191.1 -7.838 461.3 457.3 0.0 -78.36 251.5 5422.3 -46.39 213.2 -18.51 664.6 301.1 137.8 -154.3 47.67 134.8 95.18 155.11 140.9 -8.538 170.1 -20.11 -149.5 -202.3 -156.57 128.8 NaN 240.2 -273.9 254.8 -172.51 417.0 1338.0 NaN NaN NaN NaN NaN NaN NaN NaN; -30.48 1.163 -44.85 296.4 -242.8 -481.7 -330.48 -870.8 NaN NaN NaN NaN 222.1 0.0 -107.2 -41.11 -200.7 NaN 358.9 NaN -82.92 NaN NaN -99.81 30.05 NaN NaN NaN -70.14 NaN NaN NaN NaN NaN 874.19 NaN NaN NaN -164.0 NaN NaN -664.4 275.9 NaN NaN NaN NaN NaN NaN NaN; 65.33 -28.7 -22.31 223.0 -150.0 -370.3 -448.2 NaN 394.6 NaN 136.0 NaN -56.08 127.4 0.0 -189.2 138.54 431.49 147.1 NaN NaN NaN NaN 71.23 -18.93 NaN NaN NaN NaN NaN 939.07 NaN NaN NaN NaN NaN NaN 570.9 -255.22 NaN -38.77 448.1 -1327.0 NaN NaN NaN NaN NaN NaN NaN; -83.98 -25.38 -223.9 109.9 28.6 -406.8 -598.8 NaN 225.3 NaN 2889.0 NaN -194.1 38.89 865.9 0.0 287.43 NaN 1255.1 NaN -182.91 -73.85 -352.9 -262.0 -181.9 NaN NaN NaN NaN NaN NaN NaN NaN NaN 243.1 NaN NaN -196.3 22.05 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN; 1139.0 2000.0 247.5 762.8 -17.4 -118.1 -341.6 -253.1 -450.3 NaN -294.8 NaN 285.36 -15.07 64.3 -24.46 0.0 89.7 -281.6 -396.0 287.0 -111.0 NaN 882.0 617.5 NaN -139.3 NaN NaN NaN 0.1004 NaN NaN NaN NaN NaN NaN NaN -334.4 NaN -89.42 NaN NaN NaN NaN NaN NaN NaN NaN NaN; -101.6 -47.63 31.87 49.8 -132.3 -378.2 -332.9 -341.6 29.1 NaN 8.87 554.4 -156.1 NaN -207.66 NaN 117.4 0.0 -169.7 -153.7 NaN -351.6 -114.7 -205.3 -2.17 NaN 2845.0 NaN NaN NaN NaN NaN -60.78 NaN NaN NaN 160.7 -158.8 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -136.6; 24.82 -40.62 -22.97 -138.4 185.4 162.6 242.8 NaN -287.5 224.66 -266.6 99.37 38.81 -157.3 -108.5 -446.86 777.4 134.3 0.0 205.27 4.933 -152.7 -15.62 -54.86 -4.624 -0.515 NaN 230.9 0.4604 NaN 177.5 NaN -62.17 -203.0 NaN 81.57 -55.77 NaN -151.5 NaN 120.3 NaN NaN NaN NaN NaN 16.23 NaN NaN NaN; 315.3 1264.0 62.32 89.86 -151.0 339.8 -66.17 -11.0 -297.8 -165.5 -256.3 193.9 -338.5 NaN NaN NaN 493.8 -313.5 92.07 0.0 13.41 -44.7 39.63 183.4 -79.08 NaN NaN NaN NaN -208.9 NaN 228.4 -95.0 NaN -463.6 NaN -11.16 NaN -228.0 NaN -337.0 169.3 127.2 NaN NaN -322.3 NaN NaN NaN NaN; 91.46 40.25 4.68 122.9 562.2 529.0 698.2 NaN 286.3 -47.51 35.38 NaN 225.4 131.2 NaN 151.38 429.7 NaN 54.32 519.1 0.0 108.3 249.2 62.42 153.0 32.73 86.2 450.1 59.02 65.56 NaN 2.22 344.4 NaN NaN NaN -168.2 NaN 6.57 NaN 63.67 NaN NaN NaN NaN NaN NaN NaN NaN NaN; 34.01 -23.5 121.3 140.8 527.6 669.9 708.7 1633.5 82.86 190.6 -132.9 80.99 -197.7 NaN NaN -141.4 140.8 587.3 258.6 543.3 -84.53 0.0 0.0 56.33 223.1 108.9 NaN NaN NaN 149.56 NaN 177.6 315.9 NaN 215.0 NaN -91.8 NaN -160.28 NaN -96.87 NaN NaN NaN NaN NaN 361.1 NaN NaN NaN; 36.7 51.06 288.5 69.9 742.1 649.1 826.76 NaN 552.1 242.8 176.5 235.6 -20.93 NaN NaN -293.7 NaN 18.98 74.04 504.2 -157.1 0.0 0.0 -30.1 192.1 NaN NaN 116.6 NaN -64.38 NaN 86.4 168.8 NaN 363.7 NaN 111.2 NaN NaN NaN 255.8 NaN NaN -35.68 NaN NaN NaN 565.9 NaN NaN; -78.45 160.9 -4.7 134.7 856.3 709.6 1201.0 10000.0 372.0 NaN 129.5 351.9 113.9 261.1 91.13 316.9 898.2 368.5 492.0 631.0 11.8 17.97 51.9 0.0 -75.97 490.9 534.7 132.2 NaN 546.7 NaN 247.8 146.6 NaN 337.7 369.5 187.1 215.2 498.6 NaN 256.5 NaN 233.1 NaN NaN NaN 423.1 63.95 NaN 108.5; 106.8 70.32 -97.27 402.5 325.7 612.8 -274.5 622.3 518.4 NaN -171.1 383.3 -25.15 108.5 102.2 2951.0 334.9 20.18 363.5 993.4 -129.7 -8.309 -0.2266 -248.4 0.0 132.7 2213.0 NaN NaN NaN NaN NaN 593.4 NaN 1337.37 NaN NaN NaN 5143.14 309.58 -145.1 NaN NaN -209.7 NaN NaN 434.1 NaN NaN NaN; -32.69 -1.996 10.38 -97.05 261.6 252.6 417.9 NaN -142.6 NaN 129.3 NaN -94.49 NaN NaN NaN NaN NaN 0.2827 NaN 113.0 -9.639 NaN -34.68 132.9 0.0 533.2 320.2 NaN NaN 139.8 304.3 10.17 -27.7 NaN NaN 10.76 NaN -223.1 NaN 248.4 NaN NaN NaN -218.9 NaN NaN NaN NaN -4.565; 5541.0 NaN 1824.0 -127.8 561.6 511.29 360.7 815.12 -101.5 NaN NaN NaN 220.66 NaN NaN NaN 134.9 2475.0 NaN NaN 1971.0 NaN NaN 514.6 -123.1 -85.12 0.0 NaN NaN NaN NaN 2990.0 -124.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN; -52.65 16.62 21.5 40.68 609.8 914.2 1081.0 1421.0 303.7 NaN 243.8 NaN 112.4 NaN NaN NaN NaN NaN 335.7 NaN -73.09 NaN -26.06 -60.71 NaN 277.8 NaN 0.0 NaN NaN NaN 292.7 NaN NaN NaN NaN -47.37 NaN NaN NaN 469.8 NaN NaN NaN NaN NaN NaN NaN NaN NaN; -7.481 NaN 28.41 19.56 461.6 448.6 NaN NaN 160.6 NaN NaN 201.5 63.71 106.7 NaN NaN NaN NaN 161.0 NaN -27.94 NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN 31.66 NaN NaN NaN 78.92 NaN NaN NaN NaN 1004.0 NaN NaN NaN -18.27 NaN NaN; -25.31 82.64 157.3 128.8 521.6 287.0 23.48 NaN 317.5 NaN -146.3 NaN -87.31 NaN NaN NaN NaN NaN NaN 570.6 -39.46 -116.21 48.48 -133.16 NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN 262.9 NaN NaN NaN 43.37 NaN NaN NaN NaN NaN NaN NaN NaN NaN; 140.0 NaN 221.4 150.6 267.6 240.8 -137.4 838.4 135.4 NaN 152.0 NaN 9.207 NaN -213.74 NaN 192.3 NaN 169.6 NaN NaN NaN NaN NaN NaN 481.3 NaN NaN NaN NaN 0.0 NaN NaN NaN -417.2 NaN NaN NaN 302.2 NaN 347.8 NaN NaN -262.0 NaN NaN -353.5 NaN NaN NaN; 128.0 NaN 58.68 26.41 501.3 431.3 NaN NaN 138.0 245.9 21.92 NaN 476.6 NaN NaN NaN NaN NaN NaN 616.6 179.25 -40.82 21.76 48.49 NaN 64.28 2448.0 -27.45 NaN NaN NaN 0.0 6.37 NaN NaN NaN NaN NaN NaN NaN 68.55 NaN NaN NaN NaN NaN NaN NaN NaN NaN; -31.52 174.6 -154.2 1112.0 524.9 494.7 79.18 NaN -142.6 NaN 24.37 -92.26 736.4 NaN NaN NaN NaN -42.71 136.9 5256.0 -262.3 -174.5 -46.8 77.55 -185.3 125.3 4288.0 NaN NaN NaN NaN 37.1 0.0 NaN 32.9 NaN -48.33 NaN 336.25 NaN -195.1 NaN NaN NaN NaN NaN NaN NaN NaN NaN; -72.88 41.38 -101.12 614.52 68.95 967.71 NaN NaN 443.6 -55.87 -111.45 NaN 173.77 NaN NaN NaN NaN NaN 329.1 NaN NaN NaN NaN NaN NaN 174.4 NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN 2073.0 NaN -119.8 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN; 50.49 64.07 -2.504 -143.2 -25.87 695.0 -240.0 NaN 110.4 NaN 41.57 NaN -93.51 -366.51 NaN -257.2 NaN NaN NaN -180.2 NaN -215.0 -343.6 -58.43 -334.12 NaN NaN NaN 85.7 NaN 535.8 NaN -111.2 NaN 0.0 NaN NaN NaN -97.71 NaN 153.7 NaN NaN NaN NaN NaN NaN NaN NaN NaN; -165.9 573.0 -123.6 397.4 389.3 218.8 386.6 NaN 114.55 354.0 175.5 NaN NaN NaN NaN NaN NaN NaN -42.31 NaN NaN NaN NaN -85.15 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 -208.8 NaN -8.804 NaN 423.4 NaN NaN NaN NaN NaN NaN NaN NaN NaN; 47.41 124.2 395.8 419.1 738.9 528.0 NaN NaN -40.9 183.8 611.3 134.5 -217.9 NaN NaN NaN NaN 281.6 335.2 898.2 383.2 301.9 -149.8 -134.2 NaN 379.4 NaN 167.9 NaN 82.64 NaN NaN 322.42 631.5 NaN 837.2 0.0 NaN 255.0 NaN 730.8 NaN NaN NaN NaN NaN NaN 2429.0 NaN NaN; -5.132 -131.7 -237.2 -157.3 649.7 645.9 NaN NaN NaN NaN NaN NaN 167.1 NaN -198.8 116.5 NaN 159.8 NaN NaN NaN NaN NaN -124.6 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 -110.65 -117.2 NaN NaN NaN 26.35 NaN NaN NaN NaN NaN NaN; -31.95 249.0 -133.9 -240.2 64.16 172.2 -287.1 NaN 97.04 13.89 -82.12 -116.7 -158.2 49.7 10.03 -185.2 343.7 NaN 150.6 -97.77 -55.21 397.24 NaN -186.7 -374.16 223.6 NaN NaN -71.0 NaN -191.7 NaN -176.26 6.699 136.6 5.15 -137.7 50.06 0.0 -5.579 72.31 NaN NaN NaN NaN NaN NaN NaN NaN NaN; 147.3 62.4 140.6 839.83 NaN NaN NaN NaN NaN NaN NaN NaN 278.15 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 33.95 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 185.6 55.8 0.0 NaN NaN NaN NaN 111.8 NaN NaN NaN NaN NaN; 529.0 1397.0 317.6 615.8 88.63 171.0 284.4 -167.3 123.4 577.5 -234.9 65.37 -247.8 NaN 284.5 NaN -22.1 NaN -61.6 1179.0 182.2 305.4 -193.0 335.7 1107.0 -124.7 NaN 885.5 NaN -64.28 -264.3 288.1 627.7 NaN -29.34 -53.91 -198.0 NaN -28.65 NaN 0.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN; -34.36 NaN 787.9 191.6 1913.0 NaN 180.2 NaN 992.4 NaN NaN NaN 448.5 961.8 1464.0 NaN NaN NaN NaN 2450.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 -2166.0 NaN NaN NaN NaN NaN NaN NaN; 110.2 NaN 234.4 221.8 84.85 NaN NaN NaN NaN NaN NaN NaN NaN -125.2 1604.0 NaN NaN NaN NaN 2496.0 NaN NaN NaN 70.81 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 745.3 0.0 NaN NaN NaN NaN NaN NaN NaN; 13.89 -16.11 -23.88 6.214 796.9 NaN 832.2 -234.7 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -196.2 NaN 161.5 NaN NaN NaN -274.1 NaN 262.0 NaN NaN NaN NaN NaN -66.31 NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN NaN; 30.74 NaN 167.9 NaN 794.4 762.7 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 844.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -32.17 NaN NaN NaN NaN 0.0 NaN NaN NaN NaN NaN; 27.97 9.755 NaN NaN 394.8 NaN -509.3 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -70.25 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN NaN NaN; -11.92 132.4 -86.88 -19.45 517.5 NaN -205.7 NaN 156.4 NaN -3.444 NaN NaN NaN NaN NaN NaN NaN 119.2 NaN NaN -194.7 NaN 3.163 7.082 NaN NaN NaN NaN NaN 515.8 NaN NaN NaN NaN NaN NaN NaN NaN NaN 101.2 NaN NaN NaN NaN NaN 0.0 NaN NaN NaN; 39.93 543.6 NaN NaN NaN 420.0 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN -363.1 -11.3 NaN NaN NaN NaN 6.971 NaN NaN NaN NaN NaN NaN NaN 148.9 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN NaN; -23.61 161.1 142.9 274.1 -61.2 -89.24 -384.3 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN 0.0 NaN; -8.479 NaN 23.93 2.845 682.5 597.8 NaN 810.5 278.8 NaN NaN NaN NaN NaN NaN NaN NaN 221.4 NaN NaN NaN NaN NaN -79.34 NaN 176.3 NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN NaN]


#psi calculation
psi=exp(-Int_GP[n,m]/T)
return psi

end 

#finds the ln(rho) values for both the pure components and the mixture by using the interraction parameters, area fraction and Qk.
function func_lnRho(info_in,T)
    LnRHO=Array{Float64}(undef,0)

    for i in eachindex(info_in.index)
     
        ind_t1=0
        ind_t2=0
        for j in eachindex(info_in.index)
            ind_t1 += info_in.AreaFrac[j]*func_intparams(T,info_in.Group[j],info_in.Group[i])
     
            ind_t21=info_in.AreaFrac[j]*func_intparams(T,info_in.Group[i],info_in.Group[j])
            ind_t22=0
            for k in eachindex(info_in.index)
                ind_t22 += info_in.AreaFrac[k]*func_intparams(T,info_in.Group[k],info_in.Group[j])  
            end
            ind_t2 += ind_t21/ind_t22
        end

        push!(LnRHO, info_in.Qk[i]*(1-log(ind_t1)-ind_t2))
        
    end 
    #append the values onto the end of the tab
    append!(info_in.LnRho,LnRHO)
    return info_in
end 


#combines lnRhos for each group from the pure component and the mixture. does this one molecule at a time
function func_lngamV2(info_ind,info_grp)

    n=info_ind.index .∈ [info_grp.index]
    idx_grp=indexin(info_ind.index[n],info_grp.index)
    idx_ind=indexin(info_ind.index[n],info_ind.index)

    lngamR=0
    for i in eachindex(idx_ind)
        fcus_idx=(info_grp.LnRho[idx_grp[i]]-info_ind.LnRho[idx_ind[i]])*info_ind.Quan[idx_ind[i]]
        lngamR += fcus_idx
    end

    return lngamR
end

#Rk data
function f_Rk()
    # [0.901;0.6743;0.447;0.2195;1.346;1.117;1.117;0.8887;0.6606;0.5312;0.3652;1.267;1.04;0.812;1.0;1.431;0.92;0.895;1.673;1.445;0.998]

    [0.9011;0.6744;0.4469;0.2195;1.3454;1.1167;1.1173;0.8886;0.6605;0.5313;0.3652;1.2663;1.0396;0.8121;1;1.4311;0.92;0.8952;1.6724;1.4457;0.998;1.9031;1.6764;1.242;1.145;0.9183;0.6908;0.9183;1.5959;1.3692;1.1417;1.4337;1.207
    ;0.9795;1.1865;0.9597;1.06;2.9993;2.8332;2.667;1.8701;1.6434;1.3013;1.528;1.4654;1.238;1.0106;2.2564;2.0606;1.8016;2.87;2.6401;3.39;1.1562;2.0086;1.7818;1.5544;1.4199;2.057;1.877;1.651;3.168;2.4088;1.264;0.9492;1.292
    ;1.0613;2.8266;2.3144;0.791;0.6948;3.0856;2.6322;1.406;1.0105;0.615;1.38;1.6035;1.4443;1.2853;1.047;1.4838;1.303;1.1044;3.981;3.0356;2.2287;2.406;1.6493;1.8174;1.967;2.1721;2.6243;1.4515;2.1905;1.9637;2.8589;2.6322
    ;2.4054;2.1226;1.8952;1.613;1.3863;1.1589;3.474;2.8569;2.6908;2.5247]

end

#Qk data
function f_Qk()
    # [0.848;0.54;0.228;0.0;1.176;0.867;0.988;0.676;0.485;0.4;0.12;0.968;0.66;0.348;1.2;1.432;1.4;0.68;1.488;1.18;0.948]

    [0.848;0.54;0.228;0;1.176;0.867;0.988;0.676;0.485;0.4;0.12;0.968;0.66;0.348;1.2;1.432;1.4;0.68;1.488;1.18;0.948;1.728;1.42;1.188;1.088;0.78;0.468;1.1;1.544;1.236;0.924;1.244;0.936;0.624;0.94;0.632;0.816;2.113;1.833;1.553;1.724
    ;1.416;1.224;1.532;1.264;0.952;0.724;1.988;1.684;1.448;2.41;2.184;2.91;0.844;1.868;1.56;1.248;1.104;1.65;1.676;1.368;2.484;2.248;0.992;0.832;1.088;0.784;2.472;2.052;0.724;0.524;2.736;2.12;1.38;0.92;0.46;1.2;1.2632;1.0063;0.7494
    ;0.4099;1.0621;0.7639;0.4657;3.2;2.644;1.916;2.116;1.416;1.648;1.828;2.1;2.376;1.248;1.796;1.488;2.428;2.12;1.812;1.904;1.592;1.368;1.06;0.748;2.796;2.14;1.86;1.58]
end

#group index data
function f_idxarr()
    #  [1.0 1.0 1.0 1.0 2.0 2.0 2.0 2.0 2.0 3.0 3.0 4.0 4.0 4.0 5.0 6.0 7.0 8.0 9.0 9.0 10.0]'
    [1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 9.0, 10.0,
    11.0, 11.0, 12.0, 13.0, 13.0, 13.0, 13.0, 14.0, 14.0, 14.0, 15.0, 15.0, 15.0, 16.0, 16.0, 17.0, 18.0, 18.0,
    18.0, 19.0, 19.0, 20.0, 20.0, 21.0, 21.0, 21.0, 22.0, 22.0, 22.0, 23.0, 23.0, 24.0, 25.0, 26.0, 26.0, 26.0,
    27.0, 28.0, 29.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0, 34.0, 35.0, 36.0, 37.0, 38.0, 39.0, 39.0, 40.0, 40.0,
    40.0, 41.0, 42.0, 42.0, 42.0, 42.0, 43.0, 43.0, 43.0, 44.0, 45.0, 45.0, 45.0, 45.0, 45.0, 45.0, 45.0, 45.0,
    46.0, 46.0, 46.0, 46.0, 46.0, 46.0, 47.0, 47.0, 48.0, 48.0, 48.0, 49.0, 50.0, 50.0, 50.0]

end 

#adds zeros onto M_lst entries
function func_M_lst(M_lst)
    M_lst2=zeros(length(M_lst),108)
    for i in eachindex(M_lst)
        M_fcus=vec(M_lst[i])
        append!(M_fcus,(zeros(108-length(M_fcus))))
        M_lst2[i,:]=M_fcus'
    end
    return M_lst2
end 

#gives all possible 3 component compositions
function func_molfracsplit(X_arr)
    x1=Array{Float64}(undef,0)
    x2=Array{Float64}(undef,0)
    for i in eachindex(X_arr)
        append!(x1,(ones(length(X_arr)-i+1)*X_arr[i]))
        append!(x2,X_arr[1:(length(X_arr)-(i-1))])
    end

    x3=ones(length(x1))-x1-x2
    x3=round.(x3,digits=10)

    return (x1,x2,x3)
end

#molecule library
function recipies()

    name_arr=Array{Matrix}(undef,0)

    label_alkane=["ethane","propane","butane","pentane","hexane","heptane","octane","nonane","decane"]

    #alkanes
    push!(name_arr,[2 0]) #ethane
    push!(name_arr,[2 1]) #propane
    push!(name_arr,[2 2]) #butane
    push!(name_arr,[2 3]) #butane
    push!(name_arr,[2 4]) #pentane
    push!(name_arr,[2 5]) #hexane
    push!(name_arr,[2 6]) #octane
    push!(name_arr,[2 7]) #nonane
    push!(name_arr,[2 8]) #decane

    # #alkanes

    # #alkenes
    # #Aromatics
    # name_arr[1]=[0 0 0 0 0 0 0 0 0 6 0 0 0 0 0 0 0 0 0 0 0] #benzene
    # #Ketones
    # name_arr[2]=[1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0] #acetone
    # #Aldheydes
    # #Alcohols
    # name_arr[3]=[1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0 0 0] #ethanol
    # #other
    # name_arr[4]=[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 0 0 0 0] #water
    

    return name_arr

end

end