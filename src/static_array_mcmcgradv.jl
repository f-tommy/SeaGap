#using Dates
#using DelimitedFiles
#using Printf
#using Random
#using Dierckx
#using Optim
#using Statistics
#using LinearAlgebra
#using Base.Threads

export static_array_mcmcgradv
"""
    static_array_mcmcgradv(lat,dep,TR_DEPTH,NPB1,NPB2,NPB3,NPB4; gm,gs,dm,ds,rs,gmode,dmode,rmode,kernel1,kernel2,hemode,hgmode,hsmode,hdmode,he,hg,hs,hd,zest,random,NSB,fn1,fn2,fn3,fn4,nloop,nburn,NA,fno0,fno1,fno2,fno3,fno4,fno5,fno6,fno7,fno8,fno9)

Perform MCMC-based static array positioning considering sloping sound speed structure.
Prior distributions for shallow gradients and gradient depth are applied.

* `lat`: Site latitude
* `dep`: Water depth of the site
* `TR_DEPTH`: Transducer depth from the sea-surface
* `NPB1`: Number of L-NTD B-spline bases
* `NPB2`: Number of S-NTD B-spline bases
* `NPB3`: Number of shallow gradient B-spline bases
* `NPB4`: Number of gradient depth B-spline bases
* `gm`: Mean value for shallow gradient constraint
* `gs`: Std value for shallow gradient constraint
* `dm`: Mean value for gradient depth constraint
* `ds`: Std value for gradient depth constraint
* `rs`: Std value for relative gradient depth constraint
* `gmode`: Mode selection for shallow gradient constraint (0:Unifrom, 1:Normal,2:Cauchy,3:Laplace)
* `dmode`: Mode selection for gradient depth constraint (0:Unifrom, 1: Uniform 0-`dep`,2:Normal,3:Truncated Normal,4:Cauchy,5:truncated cauchy,6:Laplace,7:Truncated Laplace)
* `rmode`: Mode selection for relative gradient depth constraint (1:Unifrom, 2:Normal,3:Cauchy,4:Laplace)
* `kernel1`: If 1, Damping + Laplace smoothing for S-NTD. If 0, Damping for S-NTD
* `kernel2`: If 1, Damping + Laplace smoothing for Gradient depth. If 0, Damping for Gradient depth
* `hemode`: Prior for hyper-parameter of obervational error (0:Half-Uniform,2:Half-Normal,3:Half-Cauchy)
* `hgmode`: Prior for hyper-parameter of S-Gradient (0:Half-Uniform,2:Half-Normal,3:Half-Cauchy)
* `hdmode`: Prior for hyper-parameter of gradient depth (0:Half-Uniform,2:Half-Normal,3:Half-Cauchy)
* `hsmode`: Prior for hyper-parameter of S-NTD (0:Half-Uniform,2:Half-Normal,3:Half-Cauchy)
* `he`: Std value of the Prior for hyper-parameter of obervational error as 10^`he`
* `hg`: Std value of the Prior for hyper-parameter of S-Gradient as 10^`hg`
* `hd`: Std value of the Prior for hyper-parameter of gradient depth as 10^`hd`
* `hs`: Std value of the Prior for hyper-parameter of S-NTD as 10^`hs`
* `random`: If true, initial values for S-NTD is provided randomly
* `NSB`: Number of the perturbated bases of S-NTDs for each iteration (`NSB=100` by default)
* `nloop`: Total number of the MCMC iterations (`nloop=1200000` by default)
* `nburn`: Burn-in period of the MCMC iterations (samples less than `nburn` is excluded from the final results; `nburn=200000` by default)
* `NA`: Number of the sampling interval (`NA=5` by default; if you set (`nloop=1200000`, `nburn=200000`, and `NA=5`), you can obtain (1200000-200000)/5 samples)
* `fn1`: Input file name for the ATD offset (`fn1="tr-ant.inp"` by default)
* `fn2`: Input file name for the initial seafloor transponder positions [m] (`fn2="pxp-ini.inp"` by default)
* `fn3`: Input file name for the initial sound speed profile (`fn3="ss_prof.inp"` by default)
* `fn4`: Input file name for the observation data  (`fn4="obsdata.inp"` by default)
* `fn5`: Input file name for the initial parameters  (`fn5="initial.inp"` by default)
* `fno0`: Output file name for logging  (`fno0=log.txt` by default)
* `fno1`: Output file name for the sampled parameters after the burn-in period (`fno1=sample.out` by default)
* `fno2`: Output file name for the posterior PDFs during the mcmc iteration (`fno2=mcmc.out` by default)
* `fno3`: Output file name for the mean array displacements (`fno3=position.out` by default)
* `fno4`: Output file name for the statistical values for the sampled parameters (`fno4=statistics.out` by default)
* `fno5`: Output file name for the acceptance ratios (`fno5=acceptance.out` by default)
* `fno6`: Output file name for the travel-time residuals (`fno6=residual_sdls.out` by default)
* `fno7`: Output file name for the estimated B-spline bases (`fno7=bspline_gradv.out` by default)
* `fno8`: Output file name for the gradients (`fno8=gradient.out` by default)
* `fno9`: Output file name for the reconstructed initial file (`fno9=initial.out` by default)

# Example
    static_array_mcmc_gradv(38.0,2.7,3.0,5,100,5,1,dmode=1,gmode=1,ds=0.4,NSB=50)
"""
function static_array_mcmcgradv(lat,dep,TR_DEPTH::Vector{Float64},NPB1=5::Int64,NPB2=100::Int64,NPB3=5::Int64,NPB4=1::Int64; gm=0.0,gs=5.0e-5,dm=0.65,ds=0.15,rm=0.0,rs=0.1,gmode=1::Int64,dmode=3::Int64,rmode=1::Int64,NSB=100::Int64,nloop=600000::Int64,nburn=100000::Int64,NA=5::Int64,fno0="log.txt"::String,fno1="sample.out"::String,fno2="mcmc.out"::String,fn1="tr-ant.inp"::String,fn2="pxp-ini.inp"::String,fn3="ss_prof.inp"::String,fn4="obsdata.inp"::String,fn5="initial.inp"::String,fno3="position.out"::String,fno4="statistics.out"::String,fno5="acceptance.out"::String,fno6="residual_sdls.out"::String,fno7="bspline_gradv.out"::String,fno8="gradient.out"::String,fno9="initial.out"::String,gscale=10.0,kernel1=1::Int64,kernel2=1::Int64,zest=true::Bool,hgmode=2::Int64,hemode=0::Int64,hsmode=0::Int64,hdmode=2::Int64,he=-3.5,hs=-3.5,hg=-3.5,hd=0.5,spc=false)
  println(stderr," === GNSS-A positioning: static_array_mcmcgradv  ===")
  # --- Input check
  nds0 = size(TR_DEPTH)[1]
  if dep <= 0
    error(" static_array_mcmcgradv: dep must be more than 0")
  end
  if NPB1 < 4
    error(" static_array_mcmcgradv: NPB1 must be more than 4")
  end
  if NPB2 < 5
    error(" static_array_mcmcgradv: NPB2 must be more than 5")
  end
  if NPB3 < 0
    error(" static_array_mcmcgradv: NPB3 must be more than 0")
  end
  if NPB4 < 0
    error(" static_array_mcmcgradv: NPB4 must be more than 0")
  end
  if NSB > NPB2
    NSB = NPB2
  end
  ntr = size(TR_DEPTH)[1] # Add
  # --- Start log
  time1 = now()
  place = pwd()
  open(fno0,"w") do out0 
  println(out0,time1)
  println(out0,"static_array_mcmcgradv.jl at $place")
  nth = nthreads() # number of threads
  println(out0,"Number of threads: $nth")
  # --- Set parameters
  println(stderr," --- Set parameters")
  NP0 = 11; NC = 18 # Number of fixed parameters
  println(out0,"Number_of_B-spline_knots_L-NTD: $NPB1")
  println(out0,"Number_of_B-spline_knots_S-NTD: $NPB2")
  println(out0,"Number_of_B-spline_knots_S-Grad: $NPB3")
  println(out0,"Number_of_B-spline_knots_G-Depth: $NPB4")
  println(out0,"Default_latitude: $lat")
  println(out0,"Site_depth: $dep")
  for n in 1:nds0
    println(out0,"TR_DEPTH-$n: $TR_DEPTH[$n]")
  end
  TR_DEPTH0 = minimum(TR_DEPTH)
  println(out0,"Height_estimation: $zest")
  println(out0,"Number_of_MCMC_loop: $nloop")
  println(out0,"Burn_in_period: $nburn")
  println(out0,"Sampling_interval: $NA")
  println(out0,"Number_of_random_sampling_bases: $NSB")
  println(out0,"Prior_kernel_model_for_S-Grad: $kernel1")
  println(out0,"Prior_kernel_model_for_G-Depth: $kernel2")
  println(out0,"Prior_mode_for_HP_S-NTD: $hsmode $hs")
  println(out0,"Prior_mode_for_HP_S-Grad: $hgmode $he")
  println(out0,"Prior_mode_for_HP_Gradient-Depth: $hdmode $hd")
  println(out0,"Prior_mode_for_HP_observational-error: $hemode $he")
  println(out0,"Prior_mode_for_S-Grad: $gmode")
  println(out0,"Prior_mode_for_Gradient-Depth: $dmode")
  println(out0,"Prior_mode_for_Gradient-components: $rmode")
  println(out0,"Prior_param_for_S-Grad: $gm $gs")
  println(out0,"Prior_param_for_Gradient-Depth: $dm $ds")
  println(out0,"Prior_param_for_Gradient-components: $rs")
  # Make Prior 
  gdist = select_dist_gmode(gmode,gm,gs)
  ddist = select_dist_dmode(dmode,dm,ds,dep)
  rdist = select_dist_rmode(rmode,rm,rs)
  hedist = select_dist_hmode(hemode,2*he)
  hgdist = select_dist_hmode(hgmode,2*hg)
  hddist = select_dist_hmode(hdmode,2*hd)
  hsdist = select_dist_hmode(hsmode,2*hs)
  # --- Read data
  println(stderr," --- Read files")
  e = read_ant(fn1)
  numk, px, py, pz = read_pxppos(fn2)
  z, v, nz_st, numz = read_prof(fn3,TR_DEPTH0)
  num, nk, tp, t1, x1, y1, z1, h1, p1, r1, t2, x2, y2, z2, h2, p2, r2, nf, ids = read_obsdata(fn4)
  NP00, a00, da0, list0 = read_initial_nolimit(fn5)
  if z[end] < maximum(abs.(pz))
    error(" static_array_mcmc_all: maximum water depth of $fn3 must be deeper than site depth of $fn2")
  end
  nds = size(e)[2]
  println(out0,"Number_of_sea-surface-platforms: $nds")
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
  smin1, smax1, ds1, tb1 = mktbasis(NPB1,t1,t2,num)
  NPBV1, id1 = retrieveb(NPB1,tb1,ds1,t1,t2,num) 
  smin2, smax2, ds2, tb2 = mktbasis(NPB2,t1,t2,num)
  if spc == false
    NPBV2 = zeros(Int64,1)
    id2 = zeros(Int64,NPB2,1)
    NPBV2[1], id2[:,1] = retrieveb(NPB2,tb2,ds2,t1,t2,num) 
  else
    NPBV2 = zeros(Int64,ntr)
    id2 = zeros(Int64,NPB2,ntr)
    for k in 1:ntr
      t1e = t1[ids .== k]
      t2e = t2[ids .== k]
      nume = size(t1e)[1]
      NPBV2[k], id2[:,k] = retrieveb(NPB2,tb2,ds2,t1e,t2e,nume) 
    end
  end
  if NPB3 >= 4
    smin3, smax3, ds3, tb3 = mktbasis(NPB3,t1,t2,num)
    NPBV3, id3 = retrieveb(NPB3,tb3,ds3,t1,t2,num) 
  else
    NPBV3 = 0
  end
  if NPB4 >= 4
    smin4, smax4, ds4, tb4 = mktbasis(NPB4,t1,t2,num)
    NPBV4, id4 = retrieveb(NPB4,tb4,ds4,t1,t2,num) 
  else
    NPBV4 = 0
  end
  # --- Initialize
  println(stderr," --- Initialize")
  NP = NP0 + NPBV1 + sum(NPBV2) + 2*NPBV3 + 2*NPBV4
  NP1 = NP0 + NPBV1 + sum(NPBV2)
  d = zeros(num); dc = zeros(num); dr = zeros(num)
  a0_tmp = a00[1:NP0+NPBV1]; da = da0[1:NP0+NPBV1]
  if zest == false
    da[3] = 0.0
  end
  listv = String[]
  a_tmp = (10.0^(a00[9])).*randn(sum(NPBV2))
  da_tmp = zeros(sum(NPBV2)) 
  da_tmp[:] .= da0[end]
  for k in 1:ntr
    for i in 1:NPBV2[k]
      push!(listv,"S-NTD_$k-$i")
    end
  end
  a0 = vcat(a0_tmp,a_tmp)
  da=vcat(da,da_tmp)
  if NPB3 >= 4
    list1 = String[]; list2 = String[]
    for i in 1:NPBV3
      push!(list1,"S-Grad-EW_$i"); push!(list2,"S-Grad-NS_$i")
    end
    for j in 1:2
      for i in 1:NPBV3
        da_tmp = da0[j+3]/gscale
        a_tmp = 0.0
        a_tmp = 10^(a0[10])*randn(1)[1]
        if Distributions.pdf(gdist,a_tmp+a0[j+3]) <= 0
          a_tmp = 0.0
        end
        push!(a0,a_tmp); push!(da,da_tmp)
      end
    end
  end
  if NPB4 >= 4
    list3 = String[]; list4 = String[]
    for i in 1:NPBV4
      push!(list3,"G-Depth-EW_$i"); push!(list4,"G-Depth-NS_$i")
    end
    for j in 1:2
      for i in 1:NPBV4
        da_tmp = da0[j+5]/gscale
        a_tmp = 0.0
        a_tmp = 10^(a0[11])*randn(1)[1]
        if Distributions.pdf(gdist,a_tmp+a0[j+5]) <= 0
          a_tmp = 0.0
        end
        push!(a0,a_tmp); push!(da,da_tmp)
      end
    end
  end
  if NPB3 >= 4 
    if NPB4 >= 4 
      list = vcat(list0[1:NP0+NPBV1],listv,list1,list2,list3,list4)
    else
      list = vcat(list0[1:NP0+NPBV1],listv,list1,list2)
    end 
  else
    if NPB4 >= 4 
      list = vcat(list0[1:NP0+NPBV1],listv,list3,list4)
    else
      list = vcat(list0[1:NP0+NPBV1],listv)
    end 
  end 
  open(fno9,"w") do out9
    for i in 1:NP
      println(out9,"$(a0[i]) $(da[i]) $(list[i])")
    end
  end
  a = copy(a0)
  # --- Earth radius is fixed
  Rg, Rl = localradius(lat)
  # --- Smoothing matrix
  Lc3 = Matrix{Float64}(I,NPBV3,NPBV3)
  Lc4 = Matrix{Float64}(I,NPBV4,NPBV4)
  if kernel1 != 0 && NPB3 >= 2
    Lg3 = smoothing_matrix_1d(NPBV3)
    Le3 = Matrix{Float64}(I,NPBV3,NPBV3)
    Lc3 = Lg3 + Le3
  end
  if kernel2 != 0 && NPB4 >= 2
    Lg4 = smoothing_matrix_1d(NPBV4)
    Le4 = Matrix{Float64}(I,NPBV4,NPBV4)
    Lc4 = Lg4 + Le4
  end
  LL3 = transpose(Lc3)*Lc3
  LL4 = transpose(Lc4)*Lc4

# --- Main Anlysis --- #
  # --- Initial calculation
  println(stderr," === Initial calculation")
  println(out0," --- Start initial calculation")
  # B-spline reproduction
  a1 = a0[NP0+1:NP0+NPBV1]
  b1, sumb1 = fill_bspline_coef(NPB1,id1,a1,0.0)
  if spc == true
    a2 = zeros(NPB2,ntr)
    b2 = zeros(NPB2,ntr)
    sumb2 = zeros(ntr)
    a2[1:NPBV2[1],1] = a0[NP0+NPBV1+1:NP0+NPBV1+NPBV2[1]]
    for k in 2:ntr
      a2[1:NPBV2[k],k] = a0[NP0+NPBV1+sum(NPBV2[1:k-1])+1:NP0+NPBV1+sum(NPBV2[1:k])]
    end
    for k in 1:ntr
      b2[:,k], sumb2[k] = fill_bspline_coef(NPB2,id2[:,k],a2[1:NPBV2[k],k],0.0)
    end
  else
    a2 = zeros(findmax(NPBV2)[1],1)
    b2 = zeros(NPB2,1)
    sumb2 = zeros(1)
    a2[:,1] = a0[NP0+NPBV1+1:NP0+NPBV1+NPBV2[1]]
    b2[:,1], sumb2[1] = fill_bspline_coef(NPB2,id2[:,1],a2[:,1],0.0)
  end
  b3 = zeros(NPB3,2)
  b4 = zeros(NPB4,2)
  if NPB3 >= 4
    for mq = 1:2
      a3 =  a0[NP1+NPBV3*(mq-1)+1:NP1+NPBV3*mq]
      b3[:,mq], sumb03 = fill_bspline_coef(NPB3,id3,a3,a0[3+mq])
    end
  end
  if NPB4 >= 4
    for mq = 1:2
      a4 =  a0[NP1+NPBV3*2+NPBV4*(mq-1)+1:NP1+NPBV3*2+NPBV4*mq]
      b4[:,mq], sumb04 = fill_bspline_coef(NPB4,id4,a4,a0[5+mq])
    end
  end
  # Initialize
  hodp = zeros(num)
  tt = zeros(num); tc1 = zeros(num); tc2 = zeros(num); vert1 = zeros(num); vert2 = zeros(num)
  hh11 = zeros(num); hh12 = zeros(num); hh21 = zeros(num); hh22 = zeros(num)
  hh1 = zeros(num); hh2 = zeros(num); tc = zeros(num); td0 = zeros(num)
  td = zeros(num); tdg1 = zeros(num); tdg2 = zeros(num)
  xd = zeros(num); yd = zeros(num); vert = zeros(num); d1 = zeros(num); d2 = zeros(num)
  gd = zeros(2); td1 = zeros(num); td2 = zeros(num); dd = zeros(num)
  td3 = zeros(num); td4 = zeros(num); td5 = zeros(num); td6 = zeros(num)
  @threads for n in 1:num
    tt[n] = (t1[n] + t2[n]) / 2.0
    td1[n] = tbspline3(tt[n],ds1,tb1,b1,NPB1)
    if spc == true
      td2[n] = tbspline3(tt[n],ds2,tb2,b2[:,ids[n]],NPB2)
    else
      td2[n] = tbspline3(tt[n],ds2,tb2,b2[:,1],NPB2)
    end
    if NPB3 >= 4
      td3[n] = tbspline3(tt[n],ds3,tb3,b3[:,1],NPB3)
      td4[n] = tbspline3(tt[n],ds3,tb3,b3[:,2],NPB3)
    else
      td3[n] = a0[4]; td4[n] = a0[5]
    end
    if NPB4 >= 4
      td5[n] = tbspline3(tt[n],ds4,tb4,b4[:,1],NPB4)
      td6[n] = tbspline3(tt[n],ds4,tb4,b4[:,2],NPB4)
    else
      td5[n] = a0[6]; td6[n] = a0[7]
    end
    # --- TT
    tc1[n], vert1[n], hh11[n], hh12[n] = xyz2ttg_rapid(px[nk[n]]+a0[1],py[nk[n]]+a0[2],pz[nk[n]]+a0[3],xd1[n],yd1[n],zd1[n],Rg,Tv0[nk[n]],Vd[nk[n]],Vr[nk[n]],tr_height,cc[nk[n],1:NC])
    tc2[n], vert2[n], hh21[n], hh22[n] = xyz2ttg_rapid(px[nk[n]]+a0[1],py[nk[n]]+a0[2],pz[nk[n]]+a0[3],xd2[n],yd2[n],zd2[n],Rg,Tv0[nk[n]],Vd[nk[n]],Vr[nk[n]],tr_height,cc[nk[n],1:NC])
    xd[n] = (xd1[n]+xd2[n])/2000 ; yd[n]=(yd1[n]+yd2[n])/2000
    vert[n] = (vert1[n] + vert2[n]) / 2.0
    hh1[n] = (hh11[n] + hh21[n]) / 2.0
    hh2[n] = (hh12[n] + hh22[n]) / 2.0
    tc[n] = tc1[n] + tc2[n]
    tdg1[n] = xd[n]*td3[n] + yd[n]*td4[n]
    tdg2[n] = hh1[n]*td5[n]/2.0*td3[n] + hh2[n]*td6[n]/2.0*td4[n]
    dd[n] = (tp[n] - tc[n])*vert[n] - td1[n] -td2[n] - tdg1[n] - tdg2[n]
    hodp[n] = dd[n]*dd[n]
  end
  hod0 = sum(hodp)
  rms0 = sqrt(hod0/num)
  hod0 = -num/2*log((10^a0[8])^2) - hod0/(2*(10^a0[8])^2)
  if spc == true
    for k in 1:ntr
      hod0 += -NPBV2[k]/2*log((10^a0[9])^2) - sumb2[k]/(2*(10^a0[9])^2)
    end
  else
    hod0 += -NPBV2[1]/2*log((10^a0[9])^2) - sumb2[1]/(2*(10^a0[9])^2)
  end
  if NPB3 >= 4
    af1 = a0[NP1+1:NP1+NPBV3]; af2 = a0[NP1+NPBV3+1:NP1+NPBV3*2]
    hod0 += -NPBV3*log(10^a0[10]) -0.5*dot(af1,LL3*af1)/((10^a0[10])^2)
    hod0 += -NPBV3*log(10^a0[10]) -0.5*dot(af2,LL3*af2)/((10^a0[10])^2)
  end
  if NPB4 >= 4
    af3 = a0[NP1+2*NPBV3+1:NP1+NPBV3*2+NPBV4]; af4 = a0[NP1+2*NPBV3+1+NPBV4:NP1+NPBV3*2+NPBV4*2]
    hod0 += -NPBV4*log(10^a0[11]) -0.5*dot(af3,LL4*af3)/((10^a0[11])^2)
    hod0 += -NPBV4*log(10^a0[11]) -0.5*dot(af4,LL4*af4)/((10^a0[11])^2)
  end
  println(stderr,"   RMS: $rms0, PDF: $hod0")
  gus = a0[6] - a0[7]
  hod0 += Distributions.logpdf(ddist,a0[6])+Distributions.logpdf(ddist,a0[7])
  hod0 += Distributions.logpdf(gdist,a0[4])+Distributions.logpdf(gdist,a0[5])
  hod0 += Distributions.logpdf(rdist,gus)
  hod0 += Distributions.logpdf(hedist,10^(2*a0[8]))
  hod0 += Distributions.logpdf(hsdist,10^(2*a0[9]))
  hod0 += Distributions.logpdf(hgdist,10^(2*a0[10]))
  hod0 += Distributions.logpdf(hddist,10^(2*a0[11]))
  println(stderr,"   RMS: $rms0, PDF: $hod0")
  println(out0,"   RMS: $rms0, PDF: $hod0")
  println(out0,"   Pos: $(a0[1:3])")
  
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
    hod = 0.0
    for n in 1:nloop
      # --- MCMC perturbation
      a = copy(a0)
      for i in 1:9
        a[i] = perturbation_param_nolimit(a0[i],da[i])
      end
      if NPB3 >= 4
         a[10] = perturbation_param_nolimit(a0[10],da[10])
      end
      if NPB4 >= 4
         a[11] = perturbation_param_nolimit(a0[11],da[11])
      end
      if NPB3 >= 4 || NPB4 >= 4
        for i in NP1+1:NP
          a[i] = perturbation_param_nolimit(a0[i],da[i])
        end
      end
      for i in NP0+1:NP0+NPBV1
        a[i] = perturbation_param_nolimit(a0[i],da[i])
      end
      nss = Vector(NP0+NPBV1+1:NP1)
      @threads for i in shuffle(nss)[1:NSB]
        a[i] = perturbation_param_nolimit(a0[i],da[i])
      end
      # --- Calculate PDF
      a1 = a[NP0+1:NP0+NPBV1]
      b1, sumb1 = fill_bspline_coef(NPB1,id1,a1,0.0)
      if spc == true
        a2 = zeros(findmax(NPBV2)[1],ntr)
        b2 = zeros(NPB2,ntr)
        sumb2 = zeros(ntr)
        a2[1:NPBV2[1],1] = a[NP0+NPBV1+1:NP0+NPBV1+NPBV2[1]]
        for k in 2:ntr
          a2[1:NPBV2[k],k] = a[NP0+NPBV1+sum(NPBV2[1:k-1])+1:NP0+NPBV1+sum(NPBV2[1:k])]
        end
        for k in 1:ntr
          b2[:,k], sumb2[k] = fill_bspline_coef(NPB2,id2[:,k],a2[1:NPBV2[k],k],0.0)
        end
      else
        a2 = zeros(findmax(NPBV2)[1],1)
        b2 = zeros(NPB2,1)
        sumb2 = zeros(1)
        a2[:,1] = a[NP0+NPBV1+1:NP0+NPBV1+NPBV2[1]]
        b2[:,1], sumb2[1] = fill_bspline_coef(NPB2,id2[:,1],a2[:,1],0.0)
      end
      b3 = zeros(NPB3,2)
      b4 = zeros(NPB4,2)
      if NPB3 >= 4
        for mq = 1:2
          a3 =  a[NP1+NPBV3*(mq-1)+1:NP1+NPBV3*mq]
          b3[:,mq], sumb03 = fill_bspline_coef(NPB3,id3,a3,a[3+mq])
        end
      end
      if NPB4 >= 4
        for mq = 1:2
          a4 =  a[NP1+NPBV3*2+NPBV4*(mq-1)+1:NP1+NPBV3*2+NPBV4*mq]
          b4[:,mq], sumb04 = fill_bspline_coef(NPB4,id4,a4,a[5+mq])
        end
      end
      hodp = zeros(num)
      @threads for i in 1:num
        tt[i] = (t1[i] + t2[i]) / 2.0
        td1[i] = tbspline3(tt[i],ds1,tb1,b1,NPB1)
        if spc == true
          td2[i] = tbspline3(tt[i],ds2,tb2,b2[:,ids[i]],NPB2)
        else
          td2[i] = tbspline3(tt[i],ds2,tb2,b2[:,1],NPB2)
        end
        if NPB3 >= 4
          td3[i] = tbspline3(tt[i],ds3,tb3,b3[:,1],NPB3)
          td4[i] = tbspline3(tt[i],ds3,tb3,b3[:,2],NPB3)
        else
          td3[i] = a[4]; td4[i] = a[5]
        end
        if NPB4 >= 4
          td5[i] = tbspline3(tt[i],ds4,tb4,b4[:,1],NPB4)
          td6[i] = tbspline3(tt[i],ds4,tb4,b4[:,2],NPB4)
        else
          td5[i] = a[6]; td6[i] = a[7]
        end
        # --- TT
        tc1[i], vert1[i], hh11[i], hh12[i] = xyz2ttg_rapid(px[nk[i]]+a[1],py[nk[i]]+a[2],pz[nk[i]]+a[3],xd1[i],yd1[i],zd1[i],Rg,Tv0[nk[i]],Vd[nk[i]],Vr[nk[i]],tr_height,cc[nk[i],1:NC])
        tc2[i], vert2[i], hh21[i], hh22[i] = xyz2ttg_rapid(px[nk[i]]+a[1],py[nk[i]]+a[2],pz[nk[i]]+a[3],xd2[i],yd2[i],zd2[i],Rg,Tv0[nk[i]],Vd[nk[i]],Vr[nk[i]],tr_height,cc[nk[i],1:NC])
        xd[i] = (xd1[i]+xd2[i])/2000 ; yd[i] = (yd1[i]+yd2[i])/2000
        vert[i] = (vert1[i] + vert2[i]) / 2.0
        hh1[i] = (hh11[i] + hh21[i]) / 2.0
        hh2[i] = (hh12[i] + hh22[i]) / 2.0
        tc[i] = tc1[i] + tc2[i]
        tdg1[i] = xd[i]*td3[i] + yd[i]*td4[i]
        tdg2[i] = hh1[i]*td5[i]/2.0*td3[i] + hh2[i]*td6[i]/2.0*td4[i]
        dd[i] = (tp[i] - tc[i])*vert[i] - td1[i] - td2[i] - tdg1[i] - tdg2[i]
        hodp[i] = dd[i]*dd[i]
      end
      hod = sum(hodp)
      rms = sqrt(hod/num)
      hod = -num*log(10^a[8]) - hod/(2*(10^a[8])^2)
      if spc == true  
        for k in 1:ntr
          hod += -NPBV2[k]*log(10^a[9]) - sumb2[k]/(2*(10^a[9])^2)
        end
      else
        hod += -NPBV2[1]*log(10^a[9]) - sumb2[1]/(2*(10^a[9])^2)
      end
      if NPB3 >= 4
        af1 = a[NP1+1:NP1+NPBV3]; af2 = a[NP1+NPBV3+1:NP1+NPBV3*2]
        hod += -NPBV3*log(10^a[10]) -0.5*dot(af1,LL3*af1)/((10^a[10])^2)
        hod += -NPBV3*log(10^a[10]) -0.5*dot(af2,LL3*af2)/((10^a[10])^2)
      end
     if NPB4 >= 4
       af3 = a[NP1+2*NPBV3+1:NP1+NPBV3*3]; af4 = a[NP1+3*NPBV3+1:NP1+NPBV3*4]
       hod += -NPBV4*log(10^a[11]) -0.5*dot(af3,LL4*af3)/((10^a[11])^2)
       hod += -NPBV4*log(10^a[11]) -0.5*dot(af4,LL4*af4)/((10^a[11])^2)
      end
      gus = a[6] - a[7]
      hod += Distributions.logpdf(ddist,a[6])+Distributions.logpdf(ddist,a[7])
      hod += Distributions.logpdf(gdist,a[4])+Distributions.logpdf(gdist,a[5])
      hod += Distributions.logpdf(rdist,gus)
      hod += Distributions.logpdf(hedist,10^(2*a[8]))
      hod += Distributions.logpdf(hsdist,10^(2*a[9]))
      hod += Distributions.logpdf(hgdist,10^(2*a[10]))
      hod += Distributions.logpdf(hddist,10^(2*a[11]))
      # --- Acceptance
      acp0 = log(rand())
      acp = hod - hod0
      # --- Save
      idef = 0
      if acp >= acp0
        idef = 1
        a0 = copy(a)
        hod0 = hod
        rms0 = rms
      end
      hod = hod0
      # --- Output
      if divrem(n,NA)[2] == 0
        println(stderr," Loop $n: $idef $hod $rms0 $(a[6]) $(a[7])")
        println(out2,"$n $idef $hod $rms")
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
  dat2 = Int.(dat02[:,2])                                                 
  # --- Processing
  println(stderr," --- Post-Processing")
  tave = mean(t1)
  ave1 = mean(dat1,dims=1)
  med1 = median(dat1,dims=1)
  std1 = std(dat1,dims=1)
  min1 = minimum(dat1,dims=1)
  max1 = maximum(dat1,dims=1)
  out = permutedims(vcat(list,ave1,std1,med1,min1,max1))
  acr = mean(dat2)
  # --- Processing for NTD
  a = ave1
  numa = length(a)
  a1 = a[NP0+1:NP0+NPBV1]
  b1, sumb1 = fill_bspline_coef(NPB1,id1,a1,0.0)
  if spc == true
    a2 = zeros(findmax(NPBV2)[1],ntr)
    b2 = zeros(NPB2,ntr)
    sumb2 = zeros(ntr)
    a2[1:NPBV2[1],1] = a[NP0+NPBV1+1:NP0+NPBV1+NPBV2[1]]
    for k in 2:ntr
      a2[1:NPBV2[k],k] = a[NP0+NPBV1+sum(NPBV2[1:k-1])+1:NP0+NPBV1+sum(NPBV2[1:k])]
    end
    for k in 1:ntr
      b2[:,k], sumb2[k] = fill_bspline_coef(NPB2,id2[:,k],a2[1:NPBV2[k],k],0.0)
    end
  else
    a2 = zeros(findmax(NPBV2)[1],1)
    b2 = zeros(NPB2,1)
    sumb2 = zeros(1)
    a2[:,1] = a[NP0+NPBV1+1:NP0+NPBV1+NPBV2[1]]
    b2[:,1], sumb2[1] = fill_bspline_coef(NPB2,id2[:,1],a2[:,1],0.0)
  end
  b3 = zeros(NPB3,2)
  b4 = zeros(NPB4,2)
  if NPB3 >= 4
    for mq = 1:2
      a3 =  a[NP1+NPBV3*(mq-1)+1:NP1+NPBV3*mq]
      b3[:,mq], sumb03 = fill_bspline_coef(NPB3,id3,a3,a[3+mq])
    end
  end
  if NPB4 >= 4
    for mq = 1:2
      a4 =  a[NP1+NPBV3*2+NPBV4*(mq-1)+1:NP1+NPBV3*2+NPBV4*mq]
      b4[:,mq], sumb04 = fill_bspline_coef(NPB4,id4,a4,a[5+mq])
    end
  end
  tt = zeros(num); dt = zeros(num); tdg1 = zeros(num); tdg2 = zeros(num); td1 = zeros(num); td2 = zeros(num)
  for n in 1:num
    k = nk[n]  # PXP number
    tt[n] = (t1[n] + t2[n]) / 2.0
    td1[n] = tbspline3(tt[n],ds1,tb1,b1,NPB1)
    if spc == true
      td2[n] = tbspline3(tt[n],ds2,tb2,b2[:,ids[n]],NPB2)
    else
      td2[n] = tbspline3(tt[n],ds2,tb2,b2[:,1],NPB2)
    end
    if NPB3 >= 4
      td3[n] = tbspline3(tt[n],ds3,tb3,b3[:,1],NPB3)
      td4[n] = tbspline3(tt[n],ds3,tb3,b3[:,2],NPB3)
    else
      td3[n] = a[4]; td4[n] = a[5]
    end
    if NPB4 >= 4
      td5[n] = tbspline3(tt[n],ds4,tb4,b4[:,1],NPB4)
      td6[n] = tbspline3(tt[n],ds4,tb4,b4[:,2],NPB4)
    else
      td5[n] = a[6]; td6[n] = a[7]
    end
    # --- Calculate TT
    tc1, vert1, hh11, hh12 = xyz2ttg_rapid(px[k]+a[1],py[k]+a[2],pz[k]+a[3],xd1[n],yd1[n],zd1[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
    tc2, vert2, hh21, hh22 = xyz2ttg_rapid(px[k]+a[1],py[k]+a[2],pz[k]+a[3],xd2[n],yd2[n],zd2[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
    xd = (xd1[n]+xd2[n])/2000 ; yd = (yd1[n]+yd2[n])/2000
    vert = (vert1 + vert2) / 2.0
    hh1 = (hh11 + hh21) / 2.0
    hh2 = (hh12 + hh22) / 2.0
    tc = tc1 + tc2
    dt[n] = (tp[n] - tc)*vert
    tdg1[n] = xd*td3[n] + yd*td4[n]
    tdg2[n] = hh1*td5[n]/2.0*td3[n] + hh2*td6[n]/2.0*td4[n]
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
    println(out5,"Acceptance_ratio: $acr")
  end
  open(fno6,"w") do out6
    Base.print_array(out6,hcat(tt,nk,(xd1+xd2)/2,(yd1+yd2)/2,(zd1+zd2)/2,dt,td1+td2+tdg1+tdg2,dt-td1-td2-tdg1-tdg2,tdg1,tdg2,td1,td2,td3,td4,td5,td6,ids))  
    println(out6,"")
  end

  open(fno7,"w") do out7
    Base.print_array(out7,hcat(tt,td1,td2,td3,td4,td5,td6))  
    println(out7,"")
  end

  open(fno8,"w") do out
    gd1 = 0.5*ave1[6]*ave1[4]; gd2 = 0.5*ave1[7]*ave1[5]
    gd = (a[6]*(std1[6])^(-2.0)+a[7]*(std1[7])^(-2.0))/((std1[6])^(-2.0)+(std1[7])^(-2.0))
    println(out,"# S-Grad(EW), S-Grad(NS) D-Grad(EW) D-Grad(NS) Grad-Depth (Ave,EW,NS)")
    println(out,"$(a[4]) $(a[5]) $gd1 $gd2 $(a[6]) $(a[7]) $gd")
  end

# --- Close process --- #
  time2 = now()
  println(stderr," Start time:",time1)
  println(stderr," Finish time:",time2)
  println(out0,time2)
  end
end
