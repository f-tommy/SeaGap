#using Dates
#using DelimitedFiles
#using Printf
#using Random
#using Dierckx
#using Optim
#using Statistics
#using LinearAlgebra
#using Base.Threads

export static_array_mcmcgradc
"""
    static_array_mcmcgradc(lat,TR_DEPTH,NPB,ss,gd0,sgd; fn1,fn2,fn3,fn4,nloop,nburn,ndelay,NA,lscale,aventd,daave,daind,tscale,fno0,fno1,fno2,fno3,fno4,fno5,fno6,fno7)

Perform MCMC-based static array positioning considering a sloping sound speed strucre with a fixed number of temporal B-spline bases.
Shallow gradients and gradient depth are constrained assuming prior normal distributions. 

* `lat`: Site latitude
* `TR_DEPTH`: Transducer depth from the sea-surface
* `NPB`: Number of temporal B-spline bases
* `NSB`: Number of the perturbated bases for each iteration (`NSB=100` by default)
* `ss`: Scale paramter (width of a distribution) for the prior Cauchy distribution applying to the shallow gradients with location parameter of zero (`ss=3.e-4` by default)
* `gd0`: Mean of the prior normal distribution applying to the gradient depth [km] (`gd0=0.65` by default)
* `sgd`: Standard deviation of the prior normal distribution applying to the gradient depth [km] (`sgd=0.1` by default)
* `nloop`: Total number of the MCMC iterations (`nloop=1200000` by default)
* `nburn`: Burn-in period of the MCMC iterations (samples less than `nburn` is excluded from the final results; `nburn=200000` by default)
* `NA`: Number of the sampling interval (`NA=5` by default; if you set (`nloop=1200000`, `nburn=200000`, and `NA=5`), you can obtain (1200000-200000)/5 samples)
* `lscale`: Scaling factor for the step width of the long-term NTD parameters (`lscale=1.0` by default)
* `aventd`: If `aventd`=true, M-H pertubation for short-term NTDs are separated into average_NTD and individual NTD perturbation (`aventd=false` by default)
* `daave`: Step width for average NTD when (`aventd=true`)
* `daind`: Scaling factor for step width of individual NTD when (`aventd=true`)
* `tscale`: Temporal scaling for time in the polynomial functions (time [sec] is converted into [hour]/`tscale`, `tscale=10` by default)
* `fn1`: Input file name for an offset between a GNSS antenna and a transducer on a sea-surface platform [m] (`fn1="tr-ant.inp"` by default)
* `fn2`: Input file name for the initial seafloor transponder positions [m] (`fn2="pxp-ini.inp"` by default)
* `fn3`: Input file name for the initial sound speed profile (`fn3="ss_prof.inp"` by default)
* `fn4`: Input file name for the basic observational data  (`fn4="obsdata.inp"` by default)
* `fno0`: Output file name for logging  (`fno0=log.txt` by default)
* `fno1`: Output file name for the sampled parameters after the burn-in period (`fno1=sample.out` by default)
* `fno2`: Output file name for the posterior PDFs during the mcmc iteration (`fno2=mcmc.out` by default)
* `fno3`: Output file name for the mean array displacements (`fno3=position.out` by default)
* `fno4`: Output file name for the statistical values for the sampled parameters (`fno4=statistics.out` by default)
* `fno5`: Output file name for the acceptance ratios (`fno5=acceptance.out` by default)
* `fno6`: Output file name for the travel-time residuals (`fno6=residual_grad.out` by default)
* `fno7`: Output file name for the estimated B-spline bases (`fno7=bspline.out` by default)

# Example
    static_array_mcmcgradc(lat,TR_DEPTH,NPB,sgd=0.2)
""" 
function static_array_mcmcgradc(lat,TR_DEPTH::Vector{Float64},NPB=100::Int64, ss=1.5e-4, gd0=0.65,sgd=0.1; NSB=100::Int64,nloop=1200000::Int64,nburn=200000::Int64,NA=5::Int64,lscale=1.0,fno0="log.txt"::String,fno1="sample.out"::String,fno2="mcmc.out"::String,fn1="tr-ant.inp"::String,fn2="pxp-ini.inp"::String,fn3="ss_prof.inp"::String,fn4="obsdata.inp"::String,fn5="initial.inp"::String,fno3="position.out"::String,fno4="statistics.out"::String,fno5="acceptance.out"::String,fno6="residual_grad.out"::String,fno7="bspline.out"::String,aventd=false,daave=2.e-6,daind=10.0,tscale=10.0,ndelay=1::Int64)
  println(stderr," === GNSS-A positioning: static_array_mcmcgradc  ===")
  # --- Input check
  nds0 = size(TR_DEPTH)[1]
  if NPB < 1
    error(" static_array_mcmcgradc: NPB must be more than 1")
  end
  if NSB > NPB
    NSB = NPB
  end
  # --- Start log
  time1 = now()
  place = pwd()
  open(fno0,"w") do out0 
  println(out0,time1)
  println(out0,"static_array_mcmcgradc.jl at $place")
  nth = nthreads() # number of threads
  println(out0,"Number of threads: $nth")
  # --- Set parameters
  println(stderr," --- Set parameters")
  NP0 = 13; NC = 18 # Number of fixed parameters
  println(out0,"Number_of_B-spline_knots: $NPB")
  println(out0,"Default_latitude: $lat")
  for n in 1:nds0
    println(out0,"TR_DEPTH-$n: $TR_DEPTH[$n]")
  end
  TR_DEPTH0 = minimum(TR_DEPTH)
  println(out0,"Number_of_MCMC_loop: $nloop")
  println(out0,"Burn_in_period: $nburn")
  println(out0,"Sampling_interval: $NA")
  println(out0,"Number_of_random_sampling_bases: $NSB")
  println(out0,"Step_widths_for_NTD: $daave $daind")
  println(out0,"Time_scale_for_polynomial_functions: $tscale")
  println(out0,"Prior_Gradient-Depth_Mean: $gd0")
  println(out0,"Prior_Gradient-Depth_SD: $sgd")
  println(out0,"Prior_Shallow-Gradient_SD: $ss")
  # --- Read data
  println(stderr," --- Read files")
  e = read_ant(fn1)
  numk, px, py, pz = read_pxppos(fn2)
  z, v, nz_st, numz = read_prof(fn3,TR_DEPTH0)
  num, nk, tp, t1, x1, y1, z1, h1, p1, r1, t2, x2, y2, z2, h2, p2, r2, nf, ids = read_obsdata(fn4)
  NP, a0, a1, a2, da, list = read_initial(fn5)
  if z[end] < maximum(abs.(pz))
    error(" static_array_mcmcgradc: maximum water depth of $fn3 must be deeper than site depth of $fn2")
  end
  a = copy(a0)
  nds = size(e)[2]
  println(out0,"Number_of_sea-surface-platforms: $nds")
  # --- Average ntd
  if aventd == true
    da[NP0+1:NP] = da[NP0+1:NP] ./ daind
  end
# --- Formatting --- #
  println(stderr," --- Initial formatting")
  # --- Calculate TR position
  println(stderr," --- Calculate TR positions")
  xd1 = zeros(num); xd2 = zeros(num)
  yd1 = zeros(num); yd2 = zeros(num)
  zd1 = zeros(num); zd2 = zeros(num)
  for i in 1:num
    xd1[i], yd1[i], zd1[i] = anttena2tr(x1[i],y1[i],z1[i],h1[i],p1[i],r1[i],e[:,ids[i]])
    xd2[i], yd2[i], zd2[i] = anttena2tr(x2[i],y2[i],z2[i],h2[i],p2[i],r2[i],e[:,ids[i]])
  end
  # --- Set mean tr_height & TT corection
  println(stderr," --- TT corection")
  println(out0,"Travel-time correction: $NC")
  tr_height = ( mean(zd1) + mean(zd2) ) / 2.0
  println(stderr,"     tr_height:",tr_height)
  Tv0 = zeros(numk); Vd = zeros(numk); Vr = zeros(numk); cc = zeros(numk,NC)
  for k in 1:numk
    Tv0[k], Vd[k], Vr[k], cc[k,1:NC], rms = ttcorrection(px[k],py[k],pz[k],tr_height,z,v,nz_st,numz,TR_DEPTH0,lat)
    println(stderr,"     RMS for PxP-$k: ",rms)
    println(out0,"     RMS for PxP-$k: ",rms)
  end
  # --- Set B-spline function
  println(stderr," --- NTD basis")
  smin, smax, ds, tb = mktbasis(NPB,t1,t2,num)
  NPBV, id = retrieveb(NPB,tb,ds,t1,t2,num) 
  # --- Initialize
  NP = NP0 + NPBV
  d = zeros(num); dc = zeros(num); dr = zeros(num)
  # --- Earth radius is fixed
  Rg, Rl = localradius(lat)

# --- Main Anlysis --- #
  # --- Initial calculation
  println(stderr," === Initial calculation")
  println(out0," --- Start initial calculation")
  b = zeros(NPB)
  for m in 1:NPB
    if id[m] >= 1
      b[m] = a0[NP0+id[m]]
    end
  end
  hodp1 = zeros(num)
  hodp2 = zeros(num)
  tt = zeros(num); tc1 = zeros(num); tc2 = zeros(num); vert1 = zeros(num); vert2 = zeros(num)
  hh11 = zeros(num); hh12 = zeros(num); hh21 = zeros(num); hh22 = zeros(num)
  hh1 = zeros(num); hh2 = zeros(num); tc = zeros(num); td0 = zeros(num)
  td = zeros(num); tdg1 = zeros(num); tdg2 = zeros(num)
  xd = zeros(num); yd = zeros(num); vert = zeros(num); d1 = zeros(num); d2 = zeros(num)
  @threads for n in 1:num
    tt[n] = (t1[n] + t2[n]) / 2.0
    td[n] = tbspline3(tt[n],ds,tb,b,NPB)
    # --- TT
    tc1[n], vert1[n], hh11[n], hh12[n] = xyz2ttg_rapid(px[nk[n]]+a0[1],py[nk[n]]+a0[2],pz[nk[n]]+a0[3],xd1[n],yd1[n],zd1[n],Rg,Tv0[nk[n]],Vd[nk[n]],Vr[nk[n]],tr_height,cc[nk[n],1:NC])
    tc2[n], vert2[n], hh21[n], hh22[n] = xyz2ttg_rapid(px[nk[n]]+a0[1],py[nk[n]]+a0[2],pz[nk[n]]+a0[3],xd2[n],yd2[n],zd2[n],Rg,Tv0[nk[n]],Vd[nk[n]],Vr[nk[n]],tr_height,cc[nk[n],1:NC])
    xd[n] = (xd1[n]+xd2[n])/2000 ; yd[n]=(yd1[n]+yd2[n])/2000
    vert[n] = (vert1[n] + vert2[n]) / 2.0
    hh1[n] = (hh11[n] + hh21[n]) / 2.0
    hh2[n] = (hh12[n] + hh22[n]) / 2.0
    tc[n] = tc1[n] + tc2[n]
    tdg1[n] = xd[n]*a0[4]+yd[n]*a0[5]
    tdg2[n] = (hh1[n]*a0[4] + hh2[n]*a0[5])*a0[6]/2
    tt[n] = (tt[n] - smin)/3600/tscale
    td0[n] = a0[7] + tt[n]*a0[8] + a0[9]*tt[n]^2 + a0[10]*tt[n]^3 + a0[11]*tt[n]^4
    d1[n] = (tp[n] - tc[n])*vert[n] - td[n] - tdg2[n]
    d2[n] = (tp[n] - tc[n])*vert[n] - tdg2[n] - td0[n] - tdg1[n]
    hodp1[n] = d1[n]*d1[n]
    hodp2[n] = d2[n]*d2[n]
  end
  hod01 = sum(hodp1)
  hod02 = sum(hodp2)
  rms01 = sqrt(hod01/num)
  rms02 = sqrt(hod02/num)
  dhod01 = -0.5*log(2*pi*sgd^2)-(a0[6]-gd0)^2/(2*sgd^2)
  dhod02 = -0.5*log(2*pi*ss^2)-a0[4]^2/(2*ss^2)
  dhod02 += -0.5*log(2*pi*ss^2)-a0[5]^2/(2*ss^2)
  hod01 = -0.5*num*log(2*pi*(10^a0[12])^2) - hod01/(2*(10^a0[12])^2)
  hod02 = -0.5*num*log(2*pi*(10^a0[13])^2) - hod02/(2*(10^a0[13])^2)
  println(stderr,"   RMS; $rms01, $dhod01, $dhod02, aPDF1: $hod01, aPDF2: $hod02")
  println(out0,"   RMS; $rms01, $dhod01, $dhod02, aPDF1: $hod01, aPDF2: $hod02")
  println(out0,"   Pos; $(a0[1:3])")
  
  # --- Main calculation
  println(stderr," === Start loop calculation")
  println(out0," --- Start loop calculation")
  open(fno1,"w") do out1
  nlist = length(list)
  for n in 1:nlist
    print(out1,list[n]," ")
  end
  println(out1,"")
  open(fno2,"w") do out2
    hod1 = 0.0; hod2 = 0.0
    for n in 1:nloop
      # --- L-NTD scaling
      if n == ndelay
        da[7:11] = da[7:11] ./ lscale
      end
      # --- MCMC perturbation
      iact = divrem(n,2)[2] # iact = 0 or 1
      a = copy(a0)
      if iact == 0
        nss = Vector(14:NP)
        for i in shuffle(nss)[1:NSB]
          a[i] = perturbation_param(a0[i],da[i],a1[i],a2[i])
        end
        for i in [1 2 3 6]
          a[i] = perturbation_param(a0[i],da[i],a1[i],a2[i])
        end
        if n >= ndelay
          a[12] = perturbation_param(a0[12],da[12],a1[12],a2[12])
        end
        if aventd == true
          a[14:NP] = perturbation_param.(a[14:NP],daave,a1[14:NP],a2[14:NP],dr=(rand()-0.5)*2)
        end
      else
        for i in 7:11
          a[i] = perturbation_param(a0[i],da[i],a1[i],a2[i])
        end
        for i in [4 5]
          a[i] = perturbation_param(a0[i],da[i],a1[i],a2[i])
        end
        if n >= ndelay
          a[13] = perturbation_param(a0[13],da[13],a1[13],a2[13])
        end
      end
      # --- Calculate PDF
      b = zeros(NPB)
      @threads for m in 1:NPB
        if id[m] >= 1
          b[m] = a[NP0+id[m]]
        end
      end
      hodp1 = zeros(num); hodp2 = zeros(num)
      @threads for i in 1:num
        tt[i] = (t1[i] + t2[i]) / 2.0
        td[i] = tbspline3(tt[i],ds,tb,b,NPB)
        # --- TT
        tc1[i], vert1[i], hh11[i], hh12[i] = xyz2ttg_rapid(px[nk[i]]+a[1],py[nk[i]]+a[2],pz[nk[i]]+a[3],xd1[i],yd1[i],zd1[i],Rg,Tv0[nk[i]],Vd[nk[i]],Vr[nk[i]],tr_height,cc[nk[i],1:NC])
        tc2[i], vert2[i], hh21[i], hh22[i] = xyz2ttg_rapid(px[nk[i]]+a[1],py[nk[i]]+a[2],pz[nk[i]]+a[3],xd2[i],yd2[i],zd2[i],Rg,Tv0[nk[i]],Vd[nk[i]],Vr[nk[i]],tr_height,cc[nk[i],1:NC])
        xd[i] = (xd1[i]+xd2[i])/2000 ; yd[i] = (yd1[i]+yd2[i])/2000
        vert[i] = (vert1[i] + vert2[i]) / 2.0
        hh1[i] = (hh11[i] + hh21[i]) / 2.0
        hh2[i] = (hh12[i] + hh22[i]) / 2.0
        tc[i] = tc1[i] + tc2[i]
        tdg1[i] = xd[i]*a[4]+yd[i]*a[5]
        tdg2[i] = (hh1[i]*a[4] + hh2[i]*a[5])*a[6]/2
        tt[i] = (tt[i] - smin)/3600/tscale
        td0[i] = 0.0
        td0[i] = a[7] + tt[i]*a[8] + a[9]*tt[i]^2 + a[10]*tt[i]^3 + a[11]*tt[i]^4
        d1[i] = (tp[i] - tc[i])*vert[i] - td[i] - tdg2[i]
        d2[i] = (tp[i] - tc[i])*vert[i] - tdg2[i] - td0[i] - tdg1[i]
        hodp1[i] = d1[i]*d1[i]
        hodp2[i] = d2[i]*d2[i]
      end
      hod1 = sum(hodp1)
      hod2 = sum(hodp2)
      rms1 = sqrt(hod1/num)
      rms2 = sqrt(hod2/num)
      # Constraint by prior distribution
      dhod1 = -0.5*log(2*pi*sgd^2)-(a[6]-gd0)^2/(2*sgd^2)
      dhod2 = -0.5*log(2*pi*ss^2)-a[4]^2/(2*ss^2)
      dhod2 += -0.5*log(2*pi*ss^2)-a[5]^2/(2*ss^2)
      # Likelihood
      hod1 = -num/2*log(2*pi*(10^a[12])^2) - hod1/(2*(10^a[12])^2)
      hod2 = -num/2*log(2*pi*(10^a[13])^2) - hod2/(2*(10^a[13])^2)
      # --- Acceptance
      acp0 = log(rand())
      acp = hod1 + hod2 - hod01 - hod02 + dhod1 - dhod01 + dhod2 - dhod02
      # --- Save
      idef = 0
      if acp >= acp0
        idef = 1
        a0 = copy(a)
        hod01 = hod1
        hod02 = hod2
        dhod01 = dhod1
        dhod02 = dhod2
        rms01 = rms1
        rms02 = rms2
      end
      hod = hod01 + hod02 + dhod01 + dhod02
      # --- Output
      if divrem(n,NA)[2] == 0
        println(stderr," Loop $n: $iact $idef $hod $rms01 $rms02 $(a0[6])")
        println(out2,"$n $iact $idef $hod01 $hod02 $hod $rms01 $rms02 $dhod01 $dhod02")
        if n > nburn
          aout = copy(a0)
          Base.print_array(out1,transpose(aout))
          println(out1,"")
        end
      end
    end
  end
  end

# --- Post processing
  dat1, list = DelimitedFiles.readdlm(fno1,header=true)
  dat02 = DelimitedFiles.readdlm(fno2)
  dat2 = Int.(dat02[:,2:3])                                                 
  dat02 = dat2[dat2[:,1].==0,:][:,2]
  dat12 = dat2[dat2[:,1].==1,:][:,2]
  # --- Processing
  println(stderr," --- Post-Processing")
  tave = mean(t1)
  ave1 = mean(dat1,dims=1)
  med1 = median(dat1,dims=1)
  std1 = std(dat1,dims=1)
  min1 = minimum(dat1,dims=1)
  max1 = maximum(dat1,dims=1)
  out = permutedims(vcat(list,ave1,std1,med1,min1,max1))
  acr0 = mean(dat02)
  acr1 = mean(dat12)
  # --- Processing for NTD
  a = ave1
  numa = length(a)
  b = zeros(NPB)
  for m in 1:NPB
    if id[m] >= 1
      b[m] = a[NP0+id[m]]
    else
      b[m] = 0.0
    end
  end
  tt = zeros(num); dt = zeros(num); tdg1 = zeros(num); tdg2 = zeros(num); td = zeros(num); td0 = zeros(num)
  for n in 1:num
    k = nk[n]  # PXP number
    tt[n] = (t1[n] + t2[n]) / 2.0
    # --- Calculate TT
    tc1, vert1, hh11, hh12 = xyz2ttg_rapid(px[k]+a[1],py[k]+a[2],pz[k]+a[3],xd1[n],yd1[n],zd1[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
    tc2, vert2, hh21, hh22 = xyz2ttg_rapid(px[k]+a[1],py[k]+a[2],pz[k]+a[3],xd2[n],yd2[n],zd2[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
    xd = (xd1[n]+xd2[n])/2000 ; yd = (yd1[n]+yd2[n])/2000
    vert = (vert1 + vert2) / 2.0
    hh1 = (hh11 + hh21) / 2.0
    hh2 = (hh12 + hh22) / 2.0
    tc = tc1 + tc2
    dt[n] = (tp[n] - tc)*vert
    tdg1[n] = xd*a[4]+yd*a[5]
    tdg2[n] = (hh1*a[4] + hh2*a[5])*a[6]/2
    tt0 = (tt[n] - smin)/3600/tscale
    td0[n] = a[7] + tt0*a[8] + a[9]*tt0^2 + a[10]*tt0^3 + a[11]*tt0^4
    td[n] = tbspline3((t1[n]+t2[n])/2.0,ds,tb,b,NPB)
  end
  # --- Output
  println(stderr," --- Output for post-processing")
  open(fno3,"w") do out3
    println(out3,"$tave $(ave1[1]) $(ave1[2]) $(ave1[3]) $(std1[1]) $(std1[2]) $(std1[3])")
  end
  open(fno4,"w") do out4
    println(out4,"#Parameter mean std median min max")
    Base.print_array(out4,out)
    println(out4,"")
  end
  open(fno5,"w") do out5
    println(out5,"Acceptance_ratio_MCMC-1: $acr0")
    println(out5,"Acceptance_ratio_MCMC-2: $acr1")
  end
  open(fno6,"w") do out6
    Base.print_array(out6,hcat(tt,nk,(xd1+xd2)/2,(yd1+yd2)/2,(zd1+zd2)/2,dt,td,tdg2,td+tdg2,td0,tdg1,td0+tdg1,dt-td-tdg2))
    println(out6,"")
  end
  open(fno7,"w") do out7
    for i in 1:NPB
      println(out7,"$i $(id[i]) $(b[i])")
    end
  end

# --- Close process --- #
  time2 = now()
  println(stderr," Start time:",time1)
  println(stderr," Finish time:",time2)
  println(out0,time2)
  end
end
