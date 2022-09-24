module Fugacity

    const R=8.314462

    #a function to solve cubic equations based on the general form on wikipedia (which is just a modified Cardanos formula).
    function cubic_solve(coeff)

        #equation kx^3+yx^2+zx+w=0
        a=coeff[1]
        b=coeff[2]
        c=coeff[3]
        d=coeff[4]

        Δ_0=b^2-3*a*c
        Δ_1=2*b^3-9*a*b*c+27*a^2*d

        C=((Δ_1+sqrt(complex(Δ_1^2-4*Δ_0^3)))/2)^(1/3)
        ξ=(-1+sqrt(complex(-3)))/2

        x=Array{ComplexF64}(undef,0)
        for i in 0:2
            sol=-1/(3*a)*(b+ξ^i*C+Δ_0/(ξ^i*C))
            push!(x, sol)
        end

        x_r=Array{Float64}(undef,0)
        for i in eachindex(x)
            if abs(imag(x[i])) <=1e-10
                push!(x_r,real(x[i]))
            end
        end

        # return x
        return sort(x_r)

    end

    #Ideas gas equation to find the molar volume of a pure gas
    function Ideal_Vm(P,T)
        Vm=R*T/P
        return Vm
    end

    #Van der Waals EOS for a pure component to find compressibility (Z)
    function VdW_EOS(P,T,Pc,Tc)
        a=27/64*((R*Tc)^2/Pc)
        b=1/8*(R*Tc/Pc)

        A=a*P/(R*T)^2
        B=b*P/(R*T)

        Z=cubic_solve([1,-(1+B),A,-A*B])
        return Z[end]
    end

    #Redlich Kwong EOS for a pure component to find compressibility (Z)
    function RK_EOS(P,T,Pc,Tc)
        a=0.427480*((R^2*Tc^2.5)/Pc)
        b=0.086640*(R*Tc/Pc)

        A=a*P/(R^2*T^2.5)
        B=b*P/(R*T)

        Z=cubic_solve([1,-1,(A-B-B^2),-A*B])
        return Z[end]
    end

    #Soave Redlich Kwong EOS for a pure component to find compressibility (Z)
    function SRK_EOS(P,T,Pc,Tc,ω)
        a=0.427480*((R^2*Tc^2)/Pc)
        b=0.086640*(R*Tc/Pc)

        α=(1+(0.48508+1.55171*ω-0.15613*ω^2)*(1-sqrt(T/Tc)))^2

        A=α*a*P/(R^2*T^2)
        B=b*P/(R*T)

        Z=cubic_solve([1,-1,(A-B-B^2),-A*B])
        return Z[end]
    end

    #Van der Waals calculation for a component mixture.
    function VdW_EOS_Mix(xi,P,T,Pc,Tc)

    ai,bi=VdW_aibi_calc(Pc,Tc)
    a,b=VdW_ab_calc(xi,ai,bi)


    A=a*P/(R*T)^2
    B=b*P/(R*T)

    Z=cubic_solve([1,-(1+B),A,-A*B])[end]

    ln_phi=Array{Float64}(undef,length(ai))
    fug_coeff=copy(ln_phi)

    @. ln_phi=log(R*T*Z/(R*T*Z-b*P))-log(Z)-(bi*P/(R*T*Z-b)+2*P*sqrt(a*ai)/((T*R)^2*Z))
    @. fug_coeff=exp(ln_phi)

    return fug_coeff, Z

    end 

    function VdW_aibi_calc(Pc,Tc)

        ai=Array{Float64}(undef,length(Pc))
        bi=copy(ai)

        @. ai=27/64*((R*Tc)^2/Pc)
        @. bi=1/8*(R*Tc/Pc)

        return ai,bi

    end 

    function VdW_ab_calc(xi,ai,bi)
        a=sum(xi.*sqrt.(ai))^2
        b=sum(xi.*bi)

        return a,b
    end

    #model always assumes a basic VdW mixing rules, no binary coefficients.
    function EOS_Mix_full(xi,P,T,Pc,Tc,ω,mdl,kij)

        #correcting from matrix to vector form
        xi=vec(xi)
        Pc=vec(Pc)
        Tc=vec(Tc)
        ω=vec(ω)

        #finding the ai and bi terms for each component
        ai,bi=full_aibi_calc(Pc,Tc,T,ω,mdl)
        #running basic mixing rules on both ai and bi. This yields a and b. 
        a,b=full_ab_calc(xi,ai,bi,kij)

        #generating A and B for later use
        A=a*P/(R*T)^2
        B=b*P/(R*T)

        #pre-defining some later variables
        ln_phi=Array{Float64}(undef,length(ai))
        fug_coeff=copy(ln_phi)
        
        if mdl!="VdW"
            #for non VdW models, finding the compressibility (Z) and then each component fugacity coefficient.
            if mdl=="RK"
                Z=cubic_solve([1,-1,(A-B-B^2),-A*B])[end]
                δ_1=1
                δ_2=0

            elseif mdl=="SRK"
                Z=cubic_solve([1,-1,(A-B-B^2),-A*B])[end]
                δ_1=1
                δ_2=0

            elseif mdl=="PR"
                Z=cubic_solve([1,-(1-B),(A-2*B-3*B^2),-(A*B-B^2-B^3)])[end]
                δ_1=1+sqrt(2)
                δ_2=1-sqrt(2)
            else 
                Z=cubic_solve([1,-1,(A-B-B^2),-A*B])[end]
                δ_1=1
                δ_2=0

            end

            @. ln_phi=bi/b*(Z-1)-log(Z-B)-A/(B*(δ_2-δ_1))*(2*sqrt(ai/a)-bi/b)*log((Z+B*δ_2)/(Z+B*δ_2))
               
        else
            #alternatve method for VdW because the above equation does not work here.
            Z=cubic_solve([1,-(1+B),A,-A*B])[end]
            @. ln_phi=log(R*T*Z/(R*T*Z-b*P))-log(Z)-(bi*P/(R*T*Z-b)+2*P*sqrt(a*ai)/((T*R)^2*Z))

        end

        #applying exponential to get the fugacity coefficents.
        @. fug_coeff=exp(ln_phi)

        #return component fugacity coefficents and total compressibility.
        return fug_coeff, Z
    
        end 


    function full_aibi_calc(Pc,Tc,T,ω,mdl)

        ai=Array{Float64}(undef,length(Pc))
        bi=copy(ai)
        if mdl=="VdW"
            @. ai=27/64*((R*Tc)^2/Pc)
            @. bi=1/8*(R*Tc/Pc)
        elseif mdl=="RK"
            @. ai=0.42747*((R*Tc)^2/Pc)*sqrt(Tc/T)
            @. bi=0.08664*(R*Tc/Pc)
        elseif mdl=="SRK"
            α=copy(ai)
            @. α=(1+(0.48508+1.55171*ω-0.15613*ω^2)*(1-sqrt(T/Tc)))^2
            @. ai=0.42747*((R*Tc)^2/Pc)*α
            @. bi=0.08664*(R*Tc/Pc)
        elseif mdl=="PR"
            α=copy(ai)
            @. α=(1+(0.3796+1.485*ω-0.1644*ω^2+0.01667*ω^3)*(1-sqrt(T/Tc)))^2
            @. ai=0.457235*((R*Tc)^2/Pc)*α
            @. bi=0.077796*(R*Tc/Pc)
        else
            println("Model not known, defaulted to SRK, please select from: VdW,RK,SRK and PR")
            α=copy(ai)
            @. α=(1+(0.48508+1.55171*ω-0.15613*ω^2)*(1-sqrt(T/Tc)))^2
            @. ai=0.42747*((R*Tc)^2/Pc)*α
            @. bi=0.08664*(R*Tc/Pc)
        end

        return ai,bi

    end 

    function full_ab_calc(xi,ai,bi,kij)
        if kij==0
            kij=zeros(eachindex(ai),eachindex(ai))
        end 
        # a=sum(xi.*sqrt.(ai))^2
        # b=sum(xi.*bi)

        # if sum(kij.^2)>0
            a=sum((xi.*sqrt.(ai))*transpose(xi.*sqrt.(ai)).*(ones(axes(kij))-kij))
            b=sum(xi.*bi)

        # else
        #     a=sum(xi.*sqrt.(ai))^2
        #     b=sum(xi.*bi)
        # end

        return a,b
    end



end