#using Dates
#using DelimitedFiles
#using Statistics
#using LinearAlgebra

export solve2ntd
function solve2ntd(lat,XDUCER_DEPTH=3.0,NPB=100::Int64; type="mean",fn0="sample.out"::String,fno0="log.txt"::String,fno1="ntdgrad.out"::String,fno2="bspline.out"::String,fn1="tr-ant.inp"::String,fn2="pxp-ini.xyh"::String,fn3="ss_prof.zv"::String,fn4="obsdata.inp"::String)
  println(stderr," === NTD calculation  ===")
  # --- Set parameters
  println(stderr," --- Set parameters")
  NP0 = 13; NC = 18
  println(stderr,"Convergence eps: $eps")
  println(stderr,"Number of B-spline knots: $NPB")
  println(stderr,"Default latitude: $lat")
  println(stderr,"XDUCER_DEPTH: $XDUCER_DEPTH")
  # --- Read data
  println(stderr," --- Read files")
  e = read_ant(fn1)
  numk, px, py, pz = read_pxppos(fn2)
  z, v, nz_st, numz = read_prof(fn3,XDUCER_DEPTH)
  num, nk, tp, t1, x1, y1, z1, h1, p1, r1, t2, x2, y2, z2, h2, p2, r2, nf = read_obsdata(fn4)
  dat, list = DelimitedFiles.readdlm(fn0,header=true)
  a = mean(dat,dims=1)
  numa = length(a)

# --- Formatting --- #
  println(stderr," --- Initial formatting")
  # --- Calculate TR position
  println(stderr," --- Calculate TR positions")
  xd1 = zeros(num); xd2 = zeros(num)
  yd1 = zeros(num); yd2 = zeros(num)
  zd1 = zeros(num); zd2 = zeros(num)
  for i in 1:num
    xd1[i], yd1[i], zd1[i] = anttena2tr(x1[i],y1[i],z1[i],h1[i],p1[i],r1[i],e)
    xd2[i], yd2[i], zd2[i] = anttena2tr(x2[i],y2[i],z2[i],h2[i],p2[i],r2[i],e)
  end
  # --- Set mean xducer_height & TT corection
  println(stderr," --- TT corection")
  xducer_height = ( mean(zd1) + mean(zd1) ) / 2.0
  println(stderr,"     xducer_height:",xducer_height)
  Tv0 = zeros(numk); Vd = zeros(numk); Vr = zeros(numk); cc = zeros(numk,NC)
  for k in 1:numk
    Tv0[k], Vd[k], Vr[k], cc[k,1:NC], rms = ttcorrection(px[k],py[k],pz[k],xducer_height,z,v,nz_st,numz,XDUCER_DEPTH,lat)
    println(stderr,"     RMS for PxP-$k: ",rms)
  end
  # --- Set B-spline function
  println(stderr," --- NTD basis")
  smin, smax, ds, tb = mktbasis(NPB,t1,t2,num)
  NPBV, id = retrieveb(NPB,tb,ds,t1,t2,num) 

# --- Main Anlysis --- #
  println(stderr," === Calculation")
  Rg, Rl = localradius(lat)
  # --- Initialize
  NP = NP0 + NPBV
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
    tc1, vert1, hh11, hh12 = xyz2ttg_rapid(px[k]+a[1],py[k]+a[2],pz[k]+a[3],xd1[n],yd1[n],zd1[n],Rg,Tv0[k],Vd[k],Vr[k],xducer_height,cc[k,1:NC])
    tc2, vert2, hh21, hh22 = xyz2ttg_rapid(px[k]+a[1],py[k]+a[2],pz[k]+a[3],xd2[n],yd2[n],zd2[n],Rg,Tv0[k],Vd[k],Vr[k],xducer_height,cc[k,1:NC])
    xd = (xd1[n]+xd2[n])/2000 ; yd = (yd1[n]+yd2[n])/2000
    vert = (vert1 + vert2) / 2.0
    hh1 = (hh11 + hh21) / 2.0
    hh2 = (hh12 + hh22) / 2.0
    tc = tc1 + tc2
    dt[n] = (tp[n] - tc)*vert
    tdg1[n] = xd*a[4]+yd*a[5]
    tdg2[n] = (hh1*a[4] + hh2*a[5])*a[6]/2
    tt0 = (tt[n] - smin)/3600
    td0[n] = a[7] + tt0*a[8] + a[9]*tt0^2 + a[10]*tt0^3 + a[11]*tt0^4
    td[n] = tbspline3((t1[n]+t2[n])/2.0,ds,tb,b,NPB)
  end

# --- Output --- #
  open(fno1,"w") do out
    Base.print_array(out,hcat(tt,nk,dt,td,tdg2,td+tdg2,td0,tdg1,td0+tdg1,dt-td-tdg2))
    println(out,"")
  end
  open(fno2,"w") do out
    Base.print_array(out,hcat(collect(1:NPB),id,b))
    println(out,"")
  end
end
