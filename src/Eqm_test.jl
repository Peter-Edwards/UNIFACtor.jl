include("UNIFAC_func_lstV2.jl")
#include("UNIFACmod_func_listV2.jl")
include("Fugacity_func_list.jl")
include("EqmCalc_func_list.jl")

using NLopt

antoine_w=[5.1806,1723.64,-40.0739]
#antoine_e=

P=101.3e+3 #in atm
Tb_arr=EqmCalc.Antoine(P,[[4.65,1435.26,-64.848],[5.24677,1598.673,-46.424]],false)

T_k=78.2+273.15

x1=0.128
y1=0.128

x=[x1 1-x1]
y=[y1 1-y1]

#M_1=[17 1] #water
M_1=[19 1] #water (unifacmod)
M_2=[1 1;2 1;15 1]#ethanol
#M_2=[1 1;21 1]#acetone
#M_3=[1 2;2 5]#n-heptane
#M_4=[1 1;2 2;15 1]#1-propanol


M_lst=[M_1, M_2]

#B=UNIFAC.Activity(T_k, M_lst, x)
Activity_arr=UNIFAC.Activity(T_k, M_lst, vec(x))
fugacity_arr=Fugacity.EOS_Mix_full(y,P,T_k,[221.2*100e+3,63*100e+3],[647,513.9],[0.344,0.644],"SRK",0)
sat_pressure=EqmCalc.Antoine(T_k,[[5.1806,1723.64,-40.0739],[5.24677,1598.673,-46.424]],true)

pure_fug_w=Fugacity.EOS_Mix_full([1],sat_pressure[1],T_k,[221.2*100e+3],[647],[0.344],"SRK",0)
pure_fug_e=Fugacity.EOS_Mix_full([1],sat_pressure[2],T_k,[63*100e+3],[513.9],[0.644],"SRK",0)

unity_eval=(x[1]*sat_pressure[1]*Activity_arr[1]*pure_fug_w[2][1])/(P*fugacity_arr[2][1])+
   (x[2]*sat_pressure[2]*Activity_arr[2]*pure_fug_e[2][1])/(P*fugacity_arr[2][2])

unity_eval2=(y[1]*P*fugacity_arr[2][1])/(Activity_arr[1]*sat_pressure[1]*pure_fug_w[2][1])+
            (y[2]*P*fugacity_arr[2][2])/(Activity_arr[2]*sat_pressure[2]*pure_fug_e[2][1])


function bubble_pt(x,params)
    #x[1] is temp x[2] is frac of y1

    Activity_arr=UNIFAC.Activity(x[1], [[17 1],[1 1;2 1;15 1]],[params[1] (1-params[1])])
    fugacity_arr=Fugacity.EOS_Mix_full([x[2],1-x[2]],params[2],x[1],[221.2*100e+3,63*100e+3],[647,513.9],[0.344,0.644],"SRK",0)
    sat_pressure=EqmCalc.Antoine(x[1],[[5.1806,1723.64,-40.0739],[5.24677,1598.673,-46.424]],true)

    pure_fug_w=Fugacity.EOS_Mix_full([1],sat_pressure[1],x[1],[221.2*100e+3],[647],[0.344],"SRK",0)
    pure_fug_e=Fugacity.EOS_Mix_full([1],sat_pressure[2],x[1],[63*100e+3],[513.9],[0.644],"SRK",0)

    #unity=(params[1]*sat_pressure[1]*Activity_arr[1]*pure_fug_w[2][1])/(params[2]*fugacity_arr[2][1])+
    #     ((1-params[1])*sat_pressure[2]*Activity_arr[2]*pure_fug_e[2][1])/(params[2]*fugacity_arr[2][2])

    y1_eval=(params[1]*sat_pressure[1]*Activity_arr[1]*pure_fug_w[2][1])/(params[2]*fugacity_arr[2][1])
    y2_eval=((1-params[1])*sat_pressure[2]*Activity_arr[2]*pure_fug_e[2][1])/(params[2]*fugacity_arr[2][2])


    #return sqrt((1-unity)^2)
    return (x[2]-y1_eval)^2+((1-x[2])-y2_eval)^2
    #return Activity_arr,fugacity_arr,sat_pressure,pure_fug_e,pure_fug_w
end

function dew_pt(x,params)
    #x[1] is temp of x[2] is frac of x1
    y=[params[1],1-params[1]]
    P=params[2]

    Activity_arr=UNIFAC.Activity(x[1], [[17 1],[1 1;2 1;15 1]],[x[2] (1-x[2])])
    fugacity_arr=Fugacity.EOS_Mix_full([params[1],1-params[1]],params[2],x[1],[221.2*100e+3,63*100e+3],[647,513.9],[0.344,0.644],"SRK",0)
    sat_pressure=EqmCalc.Antoine(x[1],[[5.1806,1723.64,-40.0739],[5.24677,1598.673,-46.424]],true)

    pure_fug_w=Fugacity.EOS_Mix_full([1],sat_pressure[1],x[1],[221.2*100e+3],[647],[0.344],"SRK",0)
    pure_fug_e=Fugacity.EOS_Mix_full([1],sat_pressure[2],x[1],[63*100e+3],[513.9],[0.644],"SRK",0)
    #pure_fug_e=(1,[1,1])
    #pure_fug_w=(1,[1,1])
    #fugacity_arr=(1,[1,1])

    #unity2=(y[1]*P*fugacity_arr[2][1])/(Activity_arr[1]*sat_pressure[1]*pure_fug_w[2][1])+
    #       (y[2]*P*fugacity_arr[2][2])/(Activity_arr[2]*sat_pressure[2]*pure_fug_e[2][1])

    eval_x1=(y[1]*P*fugacity_arr[2][1])/(Activity_arr[1]*sat_pressure[1]*pure_fug_w[2][1])
    eval_x2=(y[2]*P*fugacity_arr[2][2])/(Activity_arr[2]*sat_pressure[2]*pure_fug_e[2][1])



    #return sqrt((1-unity2)^2)   
    return (x[2]-eval_x1)^2+((1-x[2])-eval_x2)^2  
end




function NLopt_dew(x1_val,P)
    
    params=[x1_val,P]

    function myfunc(x::Vector,grad::Vector)
        if length(grad) > 0
			      nothing
        end
        return dew_pt(x, params)
    end

    
    #opt = Opt(:GN_ORIG_DIRECT_L, 2)
    #opt = Opt(:LN_COBYLA, 2)
    #opt=Opt(:LN_NELDERMEAD,2)
    #opt=Opt(:GN_AGS, 2)
    #opt=Opt(:GN_ISRES, 2)
    opt=Opt(:LN_NEWUOA_BOUND, 2) 
    opt.lower_bounds = [345.0,0.0]
    opt.upper_bounds=[380.0,1.0]
    #opt.maxtime=0.5
    opt.ftol_abs=1e-7
    #opt.maxeval=80

    opt.min_objective = myfunc
    #inequality_constraint!(opt::Opt, iq_func, 1e-2)
    
    (minf,minx,ret) = optimize(opt, [360,0.5])
    #numevals = opt.numevals # the number of function evaluations
    #println("got $minf at $minx after $numevals iterations (returned $ret)")

return minf, minx

end


function NLopt_bub(x1_val,P)
    
    params=[x1_val,P]

    function myfunc(x::Vector,grad::Vector)
        if length(grad) > 0
			      nothing
        end
        return bubble_pt(x, params)
    end

    
    #opt = Opt(:GN_ORIG_DIRECT_L, 2)
    #opt = Opt(:LN_COBYLA, 2)
    #opt=Opt(:LN_NELDERMEAD,2)  
    #opt=Opt(:GN_AGS, 2)
    #opt=Opt(:GN_ISRES, 2)
    #opt=Opt(:LN_BOBYQA, 2) 
    #opt=Opt(:LN_SBPLX, 2)
    opt=Opt(:LN_NEWUOA_BOUND, 2) 
    opt.lower_bounds = [345.0,0.0]
    opt.upper_bounds=[380.0,1.0]
    #opt.maxtime=0.5
    opt.ftol_abs=1e-7
    #opt.maxeval=80

    opt.min_objective = myfunc
    #inequality_constraint!(opt::Opt, iq_func, 1e-2)
    
    (minf,minx,ret) = optimize(opt, [360,0.5])
    #numevals = opt.numevals # the number of function evaluations
    #println("got $minf at $minx after $numevals iterations (returned $ret)")

return minf, minx

end






n=10
T_tst_arr=LinRange(345,380,n)
y1_tst_arr=LinRange(0,1,n)
unity_arr=zeros(n,n)

function Surf_gen(T_tst_arr,y1_tst_arr,unity_arr)

for i in eachindex(T_tst_arr)
    for j in eachindex(y1_tst_arr )
        unity_arr[i,j]=dew_pt([T_tst_arr[i],y1_tst_arr[j]],[0.8,101.3e+3])
    end
end

return unity_arr

end


function gen_graph(n)
x1_arr=LinRange(0,1,n)
#bpt_arr=Array{Float64}(undef,n)
#dpt_arr=Array{Float64}(undef,n)

bpt_T=Array{Float64}(undef,n)
bpt_dp=Array{Float64}(undef,n)
for i in eachindex(x1_arr) 
    #bpt_arr[i]=NLopt_bub(x1_arr[i],101.3e+3)[2][1]-273.15
    #dpt_arr[i]=NLopt_dew(x1_arr[i],101.3e+3)[2][1]-273.15
    bpt_calc=NLopt_bub(x1_arr[i],101.3e+3)[2]
    bpt_T[i]=bpt_calc[1]-273.15
    bpt_dp[i]=bpt_calc[2]



end

#return x1_arr,bpt_arr,dpt_arr
return x1_arr,bpt_T,bpt_dp
end



function gen_graph2(n)
x1_arr=LinRange(0,1,n)
bpt_arr=Array{Float64}(undef,n)
dpt_arr=Array{Float64}(undef,n)

#bpt_T=Array{Float64}(undef,n)
#bpt_dp=Array{Float64}(undef,n)
for i in eachindex(x1_arr) 
    bpt_arr[i]=NLopt_bub(x1_arr[i],101.3e+3)[2][1]-273.15
    dpt_arr[i]=NLopt_dew(x1_arr[i],101.3e+3)[2][1]-273.15
   # bpt_calc=NLopt_bub(x1_arr[i],101.3e+3)[2]
   # bpt_T[i]=bpt_calc[1]
   # bpt_dp[i]=bpt_calc[2]



end

return x1_arr,bpt_arr,dpt_arr
#return x1_arr,bpt_T,bpt_dp
end
