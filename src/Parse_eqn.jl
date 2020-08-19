RO2_names=split(readline(joinpath(@__DIR__,"../data/RO2.csv")),",")

function stoich_str2int(stoich::String)::Int
    if stoich==""
        return 1
    else
        return parse(Int,stoich)
    end
end

function parse_reactants(file::String)#,RO2_names::Array{String,1}
    reactants_dict=Dict{String,Set{Tuple{Int,Int}}}() # Dict{Reactants,Set{Tuple{LineNum,Stoich}}}
    reactants2ind=Dict{String,Int}()
    ind2reactants=Dict{Int,String}()
    reactant_step=0
    #file="MCM_test.eqn.txt"
    num_eqns=0
    open(file,"r") do f
        #global num_eqns
        for line in eachline(f)
            ind = parse(Int,match(r"\{([0-9]+)\.\}",line,1).captures[1])
            num_eqns=max(ind,num_eqns)
            line_parts=split(strip(line)[1:end-1],":")
            eqn=strip(split(line_parts[1],"}")[2])
            lhs_rhs=[strip(i) for i in split(eqn,"=")]
            lhs=[strip(i) for i in split(lhs_rhs[1],"+")]
            rhs=[strip(i) for i in split(lhs_rhs[2],"+")]
            reac_regexp=r"([0-9]*)([A-Z][0-9A-Z]*)"
            if lhs[1]==""#in some eqns, O2 and H2O is ignored, like O+O3=
                reacs=[]
            else
                reacs=[Tuple{String,String}(match(reac_regexp,i).captures) for i in lhs]
            end
            if rhs[1]==""
                prods=[]
            else
                prods=[Tuple{String,String}(match(reac_regexp,i).captures) for i in rhs]
            end
            for (stoich,reac) in reacs
                num_stoich=stoich_str2int(stoich)
                if !haskey(reactants_dict,reac)
                    reactants_dict[reac]=Set()
                    reactant_step+=1
                    reactants2ind[reac]=reactant_step
                    ind2reactants[reactant_step]=reac
                end
                push!(reactants_dict[reac],(ind,num_stoich))
            end
            for (stoich,reac) in prods
                num_stoich=stoich_str2int(stoich)*(-1)# because it is product
                if !haskey(reactants_dict,reac)
                    reactants_dict[reac]=Set()
                    reactant_step+=1
                    reactants2ind[reac]=reactant_step
                    ind2reactants[reactant_step]=reac
                end
                push!(reactants_dict[reac],(ind,num_stoich))
            end
        end
    end

    num_reactants=length(reactants_dict)
    @assert(num_reactants==reactant_step)
    #reactants_list=[ for i in 1:num_reactants]
    reactants_inds=1:num_reactants
    #reactants2ind=Dict(reactants_list[i]=>i for i in reactants_inds)

    RO2_inds=[reactants2ind[i] for i in RO2_names if haskey(reactants_dict,i)]

    stoich_mtx=spzeros(num_reactants,num_eqns)
    reactants_mtx=spzeros(num_reactants,num_eqns)
    for i in reactants_inds
        for (eqn_ind,stoich) in reactants_dict[ind2reactants[i]]
            stoich_mtx[i,eqn_ind]+=stoich#FOR DUPLICATED REACTANTS A+A->B+B
            if stoich>0
                reactants_mtx[i,eqn_ind]+=stoich #FOR CATALYSE A+B->A+C
            end
        end
    end
    dropzeros!(stoich_mtx)
    return (stoich_mtx,reactants_mtx,RO2_inds,num_eqns,num_reactants,reactants2ind)
end

function gen_evaluate_rates(config::JlBoxConfig)
    file = config.file
    photolysis_config = config.photolysis_config
    if photolysis_config isa DiurnalPhotolysisConfig
        zenith_expr = quote
                # solar declination angle  
                dec = $(photolysis_config.declination)
                # latitude 
                lat = $(photolysis_config.latitude)
                #pi = 4.0*atan(1.0) 
                # local hour angle - representing time of day 
                lha = (1.0+ttime/4.32E+4)*pi 
                radian = 180.0/pi 
                lat = lat/radian 
                dec = dec/radian 
                #theta = acos(cos(lha)*cos(dec)*cos(lat)+sin(dec)*sin(lat)) 
                sinld = sin(lat)*sin(dec) 
                cosld = cos(lat)*cos(dec) 
                cosx = (cos(lha)*cosld)+sinld
                if cosx < 0
                    cosx = 0.#night? otherwise J would yield complex number
                end
        end
    elseif photolysis_config isa FixedPhotolysisConfig
        zenith_expr = quote cosx = $(photolysis_config.cos_zenith) end
    else
        println(config.io, "Error in Photolysis Config")
        exit(-1)
    end
    rate_expr=quote end
    open(file,"r") do f
        for line in eachline(f)
            line_parts=split(strip(line)[1:end-1],":")
            spd=strip(line_parts[2])
            #replacing 123.456D-12 to 123.456E-12
            num_regexp=r"([0-9]*\.?[0-9]+)D([\+\-]?[0-9]+)"
            spd_exprs=replace(spd,num_regexp=>s"\1E\2")
            spd_exprs=replace(spd_exprs,"TEMP"=>"temp")
            spd_exprs=replace(spd_exprs,"EXP"=>"exp")
            spd_exprs=replace(spd_exprs,"**"=>"^")
            spd_exprs=replace(spd_exprs,r"J\(([0-9]+)\)"=>s"J[\1]")
            spd_expr=Meta.parse(spd_exprs)
            ind = parse(Int,match(r"\{([0-9]+)\.\}",line,1).captures[1])
            push!(rate_expr.args,:(rate_values[$ind]=$spd_expr))
        end
    end
    return quote


function evaluate_rates!(ttime::FT,RO2::FT,H2O::FT,temp::FT,rate_values::Array{FT,1},J::Array{FT,1}) where {FT}
    # ttime: seconds
    # RO2: sum of RO2s
    # Now cycle through all of the mcm_constants_dict values.
    # Please note this includes photolysis rates that will change with time of day or condition in the chamber. You will
    # need to modify this potentially.
    # Creating reference to constant values used in rate expressions 

    # constants used in calculation of reaction rates 
    M  = 2.55E+19  #Check this against pressure assumed in model 
    N2 = 0.79*M 
    O2 = 0.2095*M 

    # kro2no : ro2      + no      = ro      + no2 
    #        : ro2      + no      = rono2 
    # iupac 1992 
    KRONO2    = 2.70E-12*exp(360.0/temp) 
    
    KRO2NO = 2.7E-12*exp(360.0/temp)  

    # kro2ho2: ro2      + ho2     = rooh    + o2 
    # mcm protocol [1997] 
    KRO2HO2   = 2.91E-13*exp(1300.0/temp) 

    # kapho2 : rcoo2    + ho2     = products 
    # mcm protocol [1997] 
    KAPHO2    = 5.2E-13*exp(980.0/temp) 

    # kapno  : rcoo2    + no      = products 
    # mej [1998] 
    KAPNO = 7.5E-12*exp(290.0/temp) 

    # kro2no3: ro2      + no3     = products 
    # mcm protocol [1997] 
    KRO2NO3   = 2.3E-12 

    # kno3al : no3      + rcho    = rcoo2   + hno3 
    # mcm protocol [1997] 
    KNO3AL    = 1.4E-12*exp(-1860.0/temp) 

    # kdec   : ro                 = products 
    # mcm protocol [1997] 
    KDEC      = 1.00E+06 

    KROSEC = 2.50E-14*exp(-300.0/temp) 

    KALKOXY=3.70E-14*exp(-460.0/temp)*O2 

    KALKPXY=1.80E-14*exp(-260.0/temp)*O2 

    KROPRIM = 2.50E-14*exp(-300.0/temp) 

    KCH3O2 = 1.03E-13*exp(365.0/temp) 

    K298CH3O2 = 3.5E-13 

    K14ISOM1 = 3.00E7*exp(-5300.0/temp)

    # ------------------------------------------------------------------- 
    # complex reactions 
    # ------------------------------------------------------------------- 

    # kfpan kbpan 
    # formation and decomposition of pan 
    # iupac 2001 (mcmv3.2) 
    kc0     = 3.28E-28*M*(temp/300.0)^(-6.87) 
    kci     = 1.125E-11*(temp/300.0)^(-1.105) 
    krc     = kc0/kci 
    fcc     = 0.30 
    nc      = 0.75-(1.27*log10(fcc)) 
    fc      = 10^(log10(fcc)/(1.0+((log10(krc))/nc)^2.0)) 
    KFPAN   = (kc0*kci)*fc/(kc0+kci) 

    kd0     = 1.10E-05*M*exp(-10100.0/temp) 
    kdi     = 1.90E+17*exp(-14100.0/temp) 
    krd     = kd0/kdi 
    fcd     = 0.30 
    ncd     = 0.75-(1.27*log10(fcd)) 
    fd      = 10.0^(log10(fcd)/(1.0+((log10(krd))/ncd)^2.0)) 
    KBPAN   = (kd0*kdi)*fd/(kd0+kdi) 

    KPPN0     = 1.7E-03*M*exp(-11280.0/temp) 
    KPPNI     = 8.3E+16*exp(-13940.0/temp) 
    KRPPN     = KPPN0/KPPNI 
    FCPPN     = 0.36 
    NCPPN     = 0.75-(1.27*log10(fcd)) 
    FPPN      = 10.0^(log10(fcd)/(1.0+((log10(krd))/ncd)^2.0)) 
    KBPPN   = (KPPN0*KPPNI)*FCPPN/(KPPN0+KPPNI) 

    # kmt01  : o        + no      = no2 
    # iupac 2001 (mcmv3.2) 
    k10     = 1.00E-31*M*(temp/300.0)^(-1.6) 

    k1i     = 5.0E-11*(temp/300.0)^(0.3) 
    kr1     = k10/k1i 
    fc1     = 0.85 
    nc1     = 0.75-(1.27*log10(fc1)) 
    f1      = 10.0^(log10(fc1)/(1.0+((log10(kr1)/nc1))^2.0)) 
    KMT01   = (k10*k1i)*f1/(k10+k1i) 

    # kmt02  : o        + no2     = no3 
    # iupac 2001 (mcmv3.2) 
    k20     = 1.30E-31*M*(temp/300.0)^(-1.5) 
    k2i     = 2.30E-11*(temp/300.0)^(0.24) 
    kr2     = k20/k2i 
    fc2     = 0.6 
    nc2     = 0.75-(1.27*log10(fc2)) 
    f2      = 10.0^(log10(fc2)/(1.0+((log10(kr2)/nc2))^2.0)) 
    KMT02   = (k20*k2i)*f2/(k20+k2i) 

    # kmt03  : no2      + no3     = n2o5 
    # iupac 2006, mcmv3.2 
    k30     = 3.60E-30*M*(temp/300.0)^(-4.1) 
    k3i     = 1.90E-12*(temp/300.0)^(0.2) 
    kr3     = k30/k3i 
    fc3     = 0.35 
    nc3     = 0.75-(1.27*log10(fc3)) 
    f3      = 10.0^(log10(fc3)/(1.0+((log10(kr3)/nc3))^2.0)) 
    KMT03   = (k30*k3i)*f3/(k30+k3i) 

    # kmt04  : n2o5               = no2     + no3 
    # iupac 2006, mcmv3.2 
    k40     = 1.30E-03*M*(temp/300.0)^(-3.5)*exp(-11000.0/temp) 
    k4i     = 9.70E+14*(temp/300.0)^(0.1)*exp(-11080.0/temp) 
    kr4     = k40/k4i 
    fc4     = 0.35 
    nc4     = 0.75-(1.27*log10(fc4)) 
    f4      = 10.0^(log10(fc4)/(1+((log10(kr4)/nc4))^2.0)) 
    KMT04   = (k40*k4i)*f4/(k40+k4i) 

    # kmt05  : oh       + co(+o2) = ho2     + co2 
    # iupac 2006 
    KMT05  = 1.44E-13*(1.0 + (M/4.2E19)) 

    # kmt06  : ho2      + ho2     = h2o2    + o2 
    # water enhancement factor 
    # iupac 1992 
    KMT06  = 1.0 + (1.40E-21*exp(2200.0/temp)*H2O) 

    # kmt06  = 1.0 + (2.00E-25*exp(4670.0/temp)*h2o) 
    # S+R 2005 values 

    # kmt07  : oh       + no      = hono 

    # iupac 2006, mcmv3.2 
    k70     = 7.40E-31*M*(temp/300.0)^(-2.4) 
    k7i     = 3.30E-11*(temp/300.0)^(-0.3) 
    kr7     = k70/k7i 
    fc7     = 0.81 
    nc7     = 0.75-(1.27*log10(fc7)) 
    f7      = 10.0^(log10(fc7)/(1+((log10(kr7)/nc7))^2.0)) 
    KMT07   = (k70*k7i)*f7/(k70+k7i) 

    # kmt08  : oh       + no2     = hno3 

    # iupac 2006, mcmv3.2 
    k80     = 3.2E-30*M*(temp/300.0)^(-4.5) 
    k8i     = 3.0E-11 
    kr8     = k80/k8i 
    fc8     = 0.41 
    nc8     = 0.75-(1.27*log10(fc8)) 
    f8      = 10.0^(log10(fc8)/(1.0+((log10(kr8)/nc8))^2.0)) 
    KMT08   = (k80*k8i)*f8/(k80+k8i) 

    # kmt09  : ho2      + no2     = ho2no2 
    # iupac 1997, mcmv3.2 
    
    k90     = 1.4E-31*M*(temp/300.0)^(-3.1) 
    k9i     = 4.0E-12 
    kr9     = k90/k9i 
    fc9     = 0.4 
    nc9     = 0.75-(1.27*log10(fc9)) 
    f9      = 10.0^(log10(fc9)/(1.0+((log10(kr9)/nc9))^2.0)) 
    KMT09   = (k90*k9i)*f9/(k90+k9i) 

    # kmt10  : ho2no2             = ho2     + no2 
    # iupac 1997, mcmv3.2 

    k100     = 4.10E-05*M*exp(-10650.0/temp) 
    k10i     = 6.0E+15*exp(-11170.0/temp) 
    kr10     = k100/k10i 
    fc10     = 0.4 
    nc10     = 0.75-(1.27*log10(fc10)) 
    f10      = 10.0^(log10(fc10)/(1.0+((log10(kr10)/nc10))^2.0)) 
    KMT10    = (k100*k10i)*f10/(k100+k10i) 

    # kmt11  : oh       + hno3    = h2o     + no3 
    # iupac 2006, mcmv3.2 

    k1     = 2.40E-14*exp(460.0/temp) 
    k3     = 6.50E-34*exp(1335.0/temp) 
    k4     = 2.70E-17*exp(2199.0/temp) 
    k2     = (k3*M)/(1.0+(k3*M/k4)) 
    KMT11  = k1 + k2 

    # kmt12 iupac 2006, mcmv3.2 

    k120 = 2.5E-31*((temp/300.0)^(-2.6))*M 
    k12i = 2.0E-12 
    kr12 = k120/k12i 
    fc12 = 0.53 
    nc12 = 0.75-(1.27*log10(fc12)) 
    f12  = 10.0^(log10(fc12)/(1.0+((log10(kr12)/nc12))^2.0)) 
    KMT12    = (k120*k12i)*f12/(k120+k12i) 

    # kmt13  : ch3o2    + no2     = ch3o2no2 
    # iupac 2006 

    k130     = 2.50E-30*((temp/300.0)^(-5.5))*M 
    k13i     = 1.80E-11 
    kr13     = k130/k13i 
    fc13     = 0.36 
    nc13     = 0.75-(1.27*log10(fc13)) 
    f13      = 10.0^(log10(fc13)/(1.0+((log10(kr13)/nc13))^2.0)) 
    KMT13    = (k130*k13i)*f13/(k130+k13i) 

    # kmt14  : ch3o2no2           = ch3o2   + no2 
    # iupac 2006, mcmv3.2 

    k140     = 9.00E-05*exp(-9690.0/temp)*M 
    k14i     = 1.10E+16*exp(-10560.0/temp) 
    kr14     = k140/k14i 
    fc14     = 0.36 
    nc14     = 0.75-(1.27*log10(fc14)) 
    f14      = 10.0^(log10(fc14)/(1.0+((log10(kr14)/nc14))^2.0)) 
    KMT14    = (k140*k14i)*f14/(k140+k14i) 

    # kmt15 iupac 2006, mcmv3.2 

    k150 = 8.60E-29*((temp/300.0)^(-3.1))*M 
    k15i = 9.00E-12*((temp/300.0)^(-0.85)) 
    kr15 = k150/k15i 
    fc15 = 0.48 
    nc15 = 0.75-(1.27*log10(fc15)) 
    f15  = 10.0^(log10(fc15)/(1.0+((log10(kr15)/nc15))^2.0)) 
    KMT15 = (k150*k15i)*f15/(k150+k15i) 

    # kmt16  :  oh  +  c3h6 
    # iupac 2006 

    k160     = 8.00E-27*((temp/300.0)^(-3.5))*M 
    k16i     = 3.00E-11*((temp/300.0)^(-1.0)) 
    kr16     = k160/k16i 
    fc16     = 0.5 
    nc16     = 0.75-(1.27*log10(fc16)) 
    f16      = 10.0^(log10(fc16)/(1.0+((log10(kr16)/nc16))^2.0)) 
    KMT16    = (k160*k16i)*f16/(k160+k16i) 

    # kmt17 iupac 2006 

    k170 = 5.00E-30*((temp/300.0)^(-1.5))*M 
    k17i = 1.00E-12 
    kr17 = k170/k17i 
    fc17 = (0.17*exp(-51.0/temp))+exp(-1.0*temp/204.0) 
    nc17 = 0.75-(1.27*log10(fc17)) 
    f17  = 10.0^(log10(fc17)/(1.0+((log10(kr17)/nc17))^2.0)) 
    KMT17 = (k170*k17i)*f17/(k170+k17i) 
    
    KMT18 = 9.5E-39*O2*exp(5270.0/temp)/(1+7.5E-29*O2*exp(5610.0/temp)) 

    # ************************************************************************ 
    # define photolysis reaction rates using derwent method from mcm2box.fac 
    # ************************************************************************ 

    $(zenith_expr.args...)
    #cosx = cos(theta) 
    secx = 1.0E+0/(cosx+1.0E-30) 

    # Data taken from photolysis.txt. Calculations done in the form of: 
    # j(k) = l(k)*cosx**( mm(k))*exp(-nn(k)*secx) 
    #J=[None]*62 
    #J          L           M          N 
    J[1]=6.073E-05*cosx^(1.743)*exp(-1.0*0.474*secx) 
    J[2]=4.775E-04*cosx^(0.298)*exp(-1.0*0.080*secx) 
    J[3]=1.041E-05*cosx^(0.723)*exp(-1.0*0.279*secx) 
    J[4]=1.165E-02*cosx^(0.244)*exp(-1.0*0.267*secx) 
    J[5]=2.485E-02*cosx^(0.168)*exp(-1.0*0.108*secx) 
    J[6]=1.747E-01*cosx^(0.155)*exp(-1.0*0.125*secx) 
    J[7]=2.644E-03*cosx^(0.261)*exp(-1.0*0.288*secx) 
    J[8]=9.312E-07*cosx^(1.230)*exp(-1.0*0.307*secx) 
    J[11]=4.642E-05*cosx^(0.762)*exp(-1.0*0.353*secx) 
    J[12]=6.853E-05*cosx^(0.477)*exp(-1.0*0.323*secx) 
    J[13]=7.344E-06*cosx^(1.202)*exp(-1.0*0.417*secx) 
    J[14]=2.879E-05*cosx^(1.067)*exp(-1.0*0.358*secx) 
    J[15]=2.792E-05*cosx^(0.805)*exp(-1.0*0.338*secx) 
    J[16]=1.675E-05*cosx^(0.805)*exp(-1.0*0.338*secx) 
    J[17]=7.914E-05*cosx^(0.764)*exp(-1.0*0.364*secx) 
    J[18]=1.482E-06*cosx^(0.396)*exp(-1.0*0.298*secx) 
    J[19]=1.482E-05*cosx^(0.396)*exp(-1.0*0.298*secx) 
    J[20]=7.600E-04*cosx^(0.396)*exp(-1.0*0.298*secx) 
    J[21]=7.992E-07*cosx^(1.578)*exp(-1.0*0.271*secx) 
    J[22]=5.804E-06*cosx^(1.092)*exp(-1.0*0.377*secx) 
    J[23]=2.4246E-06*cosx^(0.395)*exp(-1.0*0.296*secx) 
    J[24]=2.424E-06*cosx^(0.395)*exp(-1.0*0.296*secx) 
    J[31]=6.845E-05*cosx^(0.130)*exp(-1.0*0.201*secx) 
    J[32]=1.032E-05*cosx^(0.130)*exp(-1.0*0.201*secx) 
    J[33]=3.802E-05*cosx^(0.644)*exp(-1.0*0.312*secx) 
    J[34]=1.537E-04*cosx^(0.170)*exp(-1.0*0.208*secx) 
    J[35]=3.326E-04*cosx^(0.148)*exp(-1.0*0.215*secx) 
    J[41]=7.649E-06*cosx^(0.682)*exp(-1.0*0.279*secx) 
    J[51]=1.588E-06*cosx^(1.154)*exp(-1.0*0.318*secx) 
    J[52]=1.907E-06*cosx^(1.244)*exp(-1.0*0.335*secx) 
    J[53]=2.485E-06*cosx^(1.196)*exp(-1.0*0.328*secx) 
    J[54]=4.095E-06*cosx^(1.111)*exp(-1.0*0.316*secx) 
    J[55]=1.135E-05*cosx^(0.974)*exp(-1.0*0.309*secx) 
    J[56]=4.365E-05*cosx^(1.089)*exp(-1.0*0.323*secx) 
    J[57]=3.363E-06*cosx^(1.296)*exp(-1.0*0.322*secx) 
    J[61]=7.537E-04*cosx^(0.499)*exp(-1.0*0.266*secx) 
    
    $rate_expr
    return rate_values
end
    end
end
