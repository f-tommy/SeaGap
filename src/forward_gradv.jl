#using Random
#using Dates
#using Statistics
#using LinearAlgebra
#using DelimitedFiles
#include("read_gnssa.jl")
#include("traveltime.jl")
# Usage: forward_test(38, 3.0, [1000,1000,-3000],txtout=false)

export forward_gradv
"""
    forward_gradv(lat,TR_DEPTH,NPB1,NPB2,NPB3,NPB4; fn1,fn2,fn3,fn4,fno0,fno,offsetlntd,siglntd)       
Forward calculation for trave-times

* `lat`: Site latitude
* `TR_DEPTH`: Transducer depth from the sea-surface
* `NPB1`: Number of L-NTD B-spline basis
* `NPB2`: Number of S-NTD B-spline basis
* `NPB3`: Number of Shallow gradients B-spline basis
* `NPB4`: Number of gradient depth B-spline basis

# Example
    forward_gradv(38, [3.0], 5, 120, 3,3)
"""
function forward_gradv(lat, TR_DEPTH::Vector{Float64}, NPB1, NPB2, NPB3, NPB4; fn1="tr-ant.inp"::String,fn4="obsdata.inp"::String,fn2="pxp-ini.inp"::String,fn3="ss_prof.inp"::String,fno="synthetic.txt"::String,fno0="synthetic_NTD.txt",fn5="initial.inp"::String,offsetlntd=8.e-4,siglntd=3.e-4)
  println(stderr," === Generate synthetic dataset  ===")
  NP0 = 11; NC = 18 # Number of fixed parameters
  # --- Input check
  nds0 = size(TR_DEPTH)[1]
  TR_DEPTH0 = minimum(TR_DEPTH)
  if NPB1 < 1
    error(" NPB1 must be more than 0")
  end
  if NPB2 < 1
    error(" NPB2 must be more than 0")
  end                                                                           
  if NPB3 < 1
    error(" NPB3 must be more than 0")
  end
  if NPB4 < 1
    error(" NPB4 must be more than 0")
  end
  # --- Read data
  println(stderr," --- Read files")
  e = read_ant(fn1)
  numk, px, py, pz = read_pxppos(fn2)
  z, v, nz_st, numz = read_prof(fn3,TR_DEPTH0)
  num, nk, tp, t1, x1, y1, z1, h1, p1, r1, t2, x2, y2, z2, h2, p2, r2, nf, ids = read_obsdata(fn4)
  NP00, a00, da00, list00 = read_initial_nolimit(fn5)
  if z[end] < maximum(abs.(pz))
    error(" Maximum water depth of $fn3 must be deeper than site depth of $fn2")
  end
  nds = size(e)[2]

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
  tr_height = ( mean(zd1) + mean(zd2) ) / 2.0
  println(stderr,"     tr_height:",tr_height)
  Tv0 = zeros(numk); Vd = zeros(numk); Vr = zeros(numk); cc = zeros(numk,NC)
  for k in 1:numk
    Tv0[k], Vd[k], Vr[k], cc[k,1:NC], rms = ttcorrection(px[k],py[k],pz[k],tr_height,z,v,nz_st,numz,TR_DEPTH0,lat)
    println(stderr,"     RMS for PxP-$k: ",rms)
  end
  # --- Set B-spline function
  println(stderr," --- NTD basis")
  smin1, smax1, ds1, tb1 = mktbasis(NPB1,t1,t2,num)
  NPBV1, id1 = retrieveb(NPB1,tb1,ds1,t1,t2,num) 
  smin2, smax2, ds2, tb2 = mktbasis(NPB2,t1,t2,num)
  NPBV2, id2 = retrieveb(NPB2,tb2,ds2,t1,t2,num) 
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
  # --- Initialize                                                              
  println(stderr," --- Initialize")
  NP = NP0 + NPBV1 + NPBV2 + 2*NPBV3 + 2*NPBV4
  NP1 = NP0 + NPBV1 + NPBV2
  a = zeros(NP)
  a[1:NP0] = a00[1:NP0]
  alntd = randn(1)*offsetlntd
  a[NP0+1:NP0+NPBV1] = randn(NPBV1)*siglntd .+ alntd[1]
  a[NP0+NPBV1+1:NP0+NPBV1+NPBV2] = randn(NPBV2)*(10.0^(a[9])) 
  if NPB3 >= 4
    a[NP1+1:NP1+NPBV3*2] = randn(NPBV3*2)*(10.0^a[10]) 
  end
  if NPB4 >= 4
    a[NP1+2*NPBV3+1:NP] = randn(NPBV4*2)*(10.0^a[11]) 
  end
  # --- Earth radius is fixed
  Rg, Rl = localradius(lat)

  # B-spline reproduction
  a1 = a[NP0+1:NP0+NPBV1]
  b1, sumb1 = fill_bspline_coef(NPB1,id1,a1,0.0)
  a2 = a[NP0+NPBV1+1:NP0+NPBV1+NPBV2]
  b2, sumb2 = fill_bspline_coef(NPB2,id2,a2,0.0)
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
      a4 =  a[NP1+2*NPBV3+NPBV4*(mq-1)+1:NP1+NPBV3*2+NPBV4*mq]
      b4[:,mq], sumb04 = fill_bspline_coef(NPB4,id4,a4,a[5+mq])
    end
  end

  # --- Calculation
  err = 10.0^(a[8])
  tt = zeros(num); tc1 = zeros(num); tc2 = zeros(num); vert1 = zeros(num); vert2 = zeros(num)
  hh11 = zeros(num); hh12 = zeros(num); hh21 = zeros(num); hh22 = zeros(num)
  hh1 = zeros(num); hh2 = zeros(num); tc = zeros(num); td0 = zeros(num)
  td = zeros(num); tdg1 = zeros(num); tdg2 = zeros(num)
  xd = zeros(num); yd = zeros(num); vert = zeros(num)
  td1 = zeros(num); td2 = zeros(num)
  td3 = zeros(num); td4 = zeros(num); td5 = zeros(num); td6 = zeros(num)
  tp = zeros(num)
  for n in 1:num
    tt[n] = (t1[n] + t2[n]) / 2.0
    td1[n] = tbspline3(tt[n],ds1,tb1,b1,NPB1)
    td2[n] = tbspline3(tt[n],ds2,tb2,b2,NPB2)
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
    # --- TT
    tc1[n], vert1[n], hh11[n], hh12[n] = xyz2ttg_rapid(px[nk[n]]+a[1],py[nk[n]]+a[2],pz[nk[n]]+a[3],xd1[n],yd1[n],zd1[n],Rg,Tv0[nk[n]],Vd[nk[n]],Vr[nk[n]],tr_height,cc[nk[n],1:NC])
    tc2[n], vert2[n], hh21[n], hh22[n] = xyz2ttg_rapid(px[nk[n]]+a[1],py[nk[n]]+a[2],pz[nk[n]]+a[3],xd2[n],yd2[n],zd2[n],Rg,Tv0[nk[n]],Vd[nk[n]],Vr[nk[n]],tr_height,cc[nk[n],1:NC])
    xd[n] = (xd1[n]+xd2[n])/2000 ; yd[n]=(yd1[n]+yd2[n])/2000
    vert[n] = (vert1[n] + vert2[n]) / 2.0
    hh1[n] = (hh11[n] + hh21[n]) / 2.0
    hh2[n] = (hh12[n] + hh22[n]) / 2.0
    tc[n] = tc1[n] + tc2[n]
    tdg1[n] = xd[n]*td3[n] + yd[n]*td4[n]
    tdg2[n] = hh1[n]*td5[n]/2.0*td3[n] + hh2[n]*td6[n]/2.0*td4[n]
    tp[n] = tc[n] + (td1[n] + td2[n] + tdg1[n] + tdg2[n])/vert[n] + err*randn(1)[1]
  end
  
  # --- Output
  open(fno0,"w") do out
    Base.print_array(out,hcat(tt,td1,td2,td3,td4,td5,td6))
    println(out,"")
  end
  open(fno,"w") do out
    for n in 1:num
      @printf(out,"%i %2.6f %10.6f %6.6f %6.6f %6.6f %3.4f %3.4f %3.4f %10.6f %6.6f %6.6f %6.6f %3.4f %3.4f %3.4f %d %d\n",nk[n],tp[n],t1[n],x1[n],y1[n],z1[n],h1[n],p1[n],r1[n],t2[n],x2[n],y2[n],z2[n],h2[n],p2[n],r2[n],nf[n],ids[n])
    end
  end
end
