#using Dates
#using Statistics
#using LinearAlgebra

export static_array_sg
"""
    static_array_sg(lat,TR_DEPTH,alpha,NPB1,NPB2,gd; fn1,fn2,fn3,fn4,eps,ITMAX,delta_pos,fno0,fno1,fno2,fno3,fno4,fno5)

Perform static array positioning considering shallow gradient with a fixed gradient depth (gd).

* `lat`: Site latitude
* `TR_DEPTH`: Transducer depth from the sea-surface
* `alpha`: Hyper-parameter for S-NTD norms
* `NPB1`: Number of L-NTD B-spline bases
* `NPB2`: Number of S-NTD B-spline bases
* `gd`: Fixed gradient depth [km]
* `eps`: Convergence threshold (`eps=1.e-4` by default)
* `IMAX`: Maximum number of iterations (`IMAX=50` by default)
* `delta_pos`: Infinitesimal amount of the array displacements to calculate the Jacobian matrix (`delta_pos=1.e-4`)
* `fn1`: Input file name for the TR-based observational data (`fn1="obsdata_tr.inp"` by default)
* `fn2`: Input file name for the initial seafloor transponder positions [m] (`fn2="pxp-ini.inp"` by default)
* `fn3`: Input file name for the initial sound speed profile (`fn3="ss_prof.inp"` by default)
* `fno0`: Output file name for logging  (`fno0=log.txt` by default)
* `fno1`: Output file name for the estimated parameters and their stds (`fno1=solve.out` by default)
* `fno2`: Output file name for the estimated array displacement (`fno2=position.out` by default)
* `fno3`: Output file name for the residuals (`fno3=residual_sdls.out` by default)
* `fno4`: Output file name for the estimated B-spline bases (`fno4=S-NTD.out` by default)
* `fno5`: Output file name for the calculated ABIC value (`fno5=ABIC.out` by default)
* `fno6`: Output file name for the gradients (`fno6=gradient.out` by default)

# Example
    static_array_sg(36.2,0.8,-0.5,5,100,0.65)
"""
function static_array_sg(lat,TR_DEPTH=3.0,alpha=-5.0,NPB1=5::Int64,NPB2=100::Int64,gd=0.65; fn1="obsdata_tr.inp"::String,fn2="pxp-ini.inp"::String,fn3="ss_prof.inp"::String,eps=1.e-4,ITMAX=50::Int64, delta_pos=1.e-4, fno0="log.txt"::String,fno1="solve.out"::String,fno2="position.out"::String,fno3="residual_sdls.out"::String,fno4="S-NTD.out"::String,fno5="ABIC.out"::String,fno6="gradient.out"::String)
  println(stderr," === GNSS-A positioning: static_array  ===")
  # --- Input check
  if TR_DEPTH < 0
    error(" static_array_sg: TR_DEPTH must be positive")
  end
  if gd < 0
    error(" static_array_sg: gd (gradient depth) must be positive")
  end
  if NPB1 < 2
    error(" static_array_sg: NPB1 must be more than 2")
  end
  if NPB2 < 2
    error(" static_array_sg: NPB2 must be more than 2")
  end
  # --- Start log
  time1 = now()
  place = pwd()
  open(fno0,"w") do out0 
  println(out0,time1)
  println(out0,"static_array_sg.jl at $place")
  # --- Set parameters
  println(stderr," --- Set parameters")
  NP0 = 5; NC = 18 # Number of fixed parameters
  dx = delta_pos; dy = delta_pos; dz = delta_pos
  println(out0,"Convergence_eps: $eps")
  println(out0,"Fixed_gradient_depth: $gd")
  println(out0,"Hyper-parameter_for_S-NTD: 10^ $alpha")
  println(out0,"Number_of_B-spline_knots1: $NPB1")
  println(out0,"Number_of_B-spline_knots2: $NPB2")
  println(out0,"Default_latitude: $lat")
  println(out0,"Maximum_iterations: $ITMAX")
  println(out0,"Delta_position: $delta_pos")
  println(out0,"TR_DEPTH: $TR_DEPTH")
  # --- Read data
  println(stderr," --- Read files")
  num, nk, tp, t1, xd1, yd1, zd1, t2, xd2, yd2, zd2, nf, ids = read_obsdata_tr(fn1)
  numk, px, py, pz = read_pxppos(fn2)
  z, v, nz_st, numz = read_prof(fn3,TR_DEPTH)
  if z[end] < maximum(abs.(pz))                                                 
    error(" static_array: maximum water depth of $fn3 must be deeper than site depth of $fn2")
  end

# --- Formatting --- #
  println(stderr," --- Initial formatting")
  # --- Set mean tr_height & TT corection
  println(stderr," --- TT corection")
  println(out0,"Travel-time correction: $NC")
  tr_height = ( mean(zd1) + mean(zd2) ) / 2.0
  println(stderr,"     Transducer_height:",tr_height)
  Tv0 = zeros(numk); Vd = zeros(numk); Vr = zeros(numk); cc = zeros(numk,NC)
  for k in 1:numk
    Tv0[k], Vd[k], Vr[k], cc[k,1:NC], rms = ttcorrection(px[k],py[k],pz[k],tr_height,z,v,nz_st,numz,TR_DEPTH,lat)
    println(stderr,"     RMS for PxP-$k: ",rms)
    println(out0,"     RMS for PxP-$k: ",rms)
  end
  # --- Set B-spline function
  println(stderr," --- NTD basis")
  smin1, smax1, ds1, tb1 = mktbasis(NPB1,t1,t2,num)
  NPBV1, id1 = retrieveb(NPB1,tb1,ds1,t1,t2,num) 
  smin2, smax2, ds2, tb2 = mktbasis(NPB2,t1,t2,num)
  NPBV2, id2 = retrieveb(NPB2,tb2,ds2,t1,t2,num) 
  # --- Initialize
  NP = NP0 + NPBV1 + NPBV2
  d = zeros(num); H = zeros(num,NP); a0 = zeros(NP); a = zeros(NP)
  H1 = zeros(num,NP); H2 = zeros(num,NP)
  dc = zeros(num); dr = zeros(num); delta = 1.e6; rms = 1.e6
  G0 = zeros(NP,NP)
  d1 = zeros(num); d2 = zeros(num); d3 = zeros(num); d4 = zeros(num)
  xd = zeros(num); yd = zeros(num)
  sigma2 = 0.0; abic = 0.0; Hinv = zeros(NP,NP)
  Rg, Rl = localradius(lat)
  for i in NP0+NPBV1+1:NP
    G0[i,i] = 10.0^alpha
  end

# --- Main Anlysis --- #
  println(stderr," === Inversion")
  println(out0,"Start iteration")
  it = 1
  while delta > eps
    if it > ITMAX
      break
    end
    println(stderr," --- Iteration: $it")
    # --- Set H-matrix
    for n in 1:num
      k = nk[n]  # PXP number
      # --- Calculate TT
      tc1, vert1, hh11, hh12 = xyz2ttg_rapid(px[k]+a0[1],py[k]+a0[2],pz[k]+a0[3],xd1[n],yd1[n],zd1[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
      tc2, vert2, hh21, hh22 = xyz2ttg_rapid(px[k]+a0[1],py[k]+a0[2],pz[k]+a0[3],xd2[n],yd2[n],zd2[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
      vert = (vert1 + vert2) / 2.0
      hh1 = (hh11 + hh21) / 2.0
      hh2 = (hh12 + hh22) / 2.0
      tc = tc1 + tc2
      d[n] = (tp[n] - tc)*vert
      # --- Differential
      tcx1, to1, vert1 = xyz2tt_rapid(px[k]+a0[1]+dx,py[k]+a0[2],pz[k]+a0[3],xd1[n],yd1[n],zd1[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
      tcx2, to2, vert2 = xyz2tt_rapid(px[k]+a0[1]+dx,py[k]+a0[2],pz[k]+a0[3],xd2[n],yd2[n],zd2[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
      tcx = tcx1 + tcx2
      tcy1, to1, vert1 = xyz2tt_rapid(px[k]+a0[1],py[k]+a0[2]+dy,pz[k]+a0[3],xd1[n],yd1[n],zd1[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
      tcy2, to2, vert2 = xyz2tt_rapid(px[k]+a0[1],py[k]+a0[2]+dy,pz[k]+a0[3],xd2[n],yd2[n],zd2[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
      tcy = tcy1 + tcy2
      tcz1, to1, vert1 = xyz2tt_rapid(px[k]+a0[1],py[k]+a0[2],pz[k]+a0[3]+dz,xd1[n],yd1[n],zd1[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
      tcz2, to2, vert2 = xyz2tt_rapid(px[k]+a0[1],py[k]+a0[2],pz[k]+a0[3]+dz,xd2[n],yd2[n],zd2[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
      tcz = tcz1 + tcz2
      # --- Fill matrix
      xd[n] = (xd1[n]+xd2[n])/2000 ; yd[n] = (yd1[n]+yd2[n])/2000
      H[n,1] = (tcx-tc)/dx*vert; H[n,2]=(tcy-tc)/dy*vert; H[n,3]=(tcz-tc)/dz*vert
      H[n,4] = xd[n]+hh1*gd/2.0; H[n,5] = yd[n]+hh2*gd/2.0
      H1[n,4] = xd[n]; H1[n,5] = yd[n]
      H2[n,4] = hh1*gd/2.0; H2[n,5] = hh2*gd/2.0
      if it == 1
        for m in 1:NPB1
          if id1[m] >= 1
            b0 = zeros(NPB1)
            b0[m] = 1.0
            H[n,NP0+id1[m]] = tbspline3((t1[n]+t2[n])/2.0,ds1,tb1,b0,NPB1)
          end
        end
        for m in 1:NPB2
          if id2[m] >= 1
            b0 = zeros(NPB2)
            b0[m] = 1.0
            H[n,NP0+NPBV1+id2[m]] = tbspline3((t1[n]+t2[n])/2.0,ds2,tb2,b0,NPB2)
          end
        end
      end
    end
    Horg = transpose(H)*H+G0
    Hinv = inv(Horg)
    a = Hinv*transpose(H)*d
    dc = H*a
    dr = d - dc
    a1 = zeros(NP); a2 = zeros(NP); a3 = zeros(NP)
    a1[4] = a[4]; a1[5] = a[5]
    a2[NP0+1:NP0+NPBV1] = a[NP0+1:NP0+NPBV1]
    a3[NP0+NPBV1+1:NP] = a[NP0+NPBV1+1:NP]
    d1 = H1*a1; d2 = H2*a1; d3 = H*a2; d4 = H*a3
    rms = std(dr)
    sa = LinearAlgebra.dot(dr,dr)
    sigma2 = sa / (num-NP)
    abic = num*log(sa) - NPBV2*log(10.0^alpha) + logdet(Horg)
    delta = std(a[1:3])
    a0[1:3] += a[1:3]
    a0[4:NP] = a[4:NP]
    println(stderr," Temporal position: $(a0[1:3]), $delta, $rms")
    println(out0,"     Iteration: $it $(a0[1]) $(a0[2]) $(a0[3]) $delta $rms")
    it += 1
  end
  println(stderr," End of loop ",it-1)
  println(stderr," --- Final position: $(a0[1:3]), $delta, $rms")
  println(out0,"End of iteration")
  cv = sqrt.(sigma2*abs.(diag(Hinv))) # Error
  tc = (t1[1] + t2[num]) / 2.0
  a = transpose(a0[1:3])
  # --- Fill NTD basis
  b1 = zeros(NPB1)
  for m in 1:NPB1
    if id1[m] >= 1
      b1[m] = a0[NP0+id1[m]]
    else
      b1[m] = 0.0
    end
  end
  b2 = zeros(NPB2)
  for m in 1:NPB2
    if id2[m] >= 1
      b2[m] = a0[NP0+NPBV1+id2[m]]
    else
      b2[m] = 0.0
    end
  end
  td = zeros(num)
  for n in 1:num
    td[n] = tbspline3((t1[n]+t2[n])/2.0,ds1,tb1,b1,NPB1)
    td[n] += tbspline3((t1[n]+t2[n])/2.0,ds2,tb2,b2,NPB2)
  end

# --- Output --- #
  open(fno1,"w") do out
    Base.print_array(out,hcat(a0,cv))
    println(out,"")
  end
  open(fno2,"w") do out
    Base.print_array(out,hcat(tc,a,transpose(cv[1:3])))
    println(out,"")
  end
  open(fno3,"w") do out
    Base.print_array(out,hcat((t1+t2)/2.0,nk,(xd1+xd2)/2.0,(yd1+yd2)/2.0,(zd1+zd2)/2.0,d,dc,dr,d1,d2,d3,d4,ids))
    println(out,"")
  end
  open(fno4,"w") do out
    Base.print_array(out,hcat(collect(1:NPB2),id2,tb2,b2))
    println(out,"")
  end
  open(fno5,"a") do out
    println(out,"$alpha $NPB1 $NPB2 $gd $abic")
  end
  dg1=a0[4]*gd/2.0; dg2=a0[5]*gd/2.0
  open(fno6,"w") do out
    println(out,"# S-Grad(EW), S-Grad(NS) D-Grad(EW) D-Grad(NS) Grad-Depth")
    println(out,"$(a0[4]) $(a0[5]) $dg1 $dg2 $gd")
  end

# --- Close process --- #
  time2 = now()
  println(stderr," Start time:",time1)
  println(stderr," Finish time:",time2)
  println(out0,time2)
  end
end
