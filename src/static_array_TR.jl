#using Dates
#using Statistics
#using LinearAlgebra

export static_array_TR
"""
    static_array_TR(lat,TR_DEPTH,NPB; fn1,fn2,fn3,fn4,eps,ITMAX,delta_pos,delta_offset,fno0,fno1,fno2,fno3,fno4,fno5)

Perform static array positioning with optimizing offset between a GNSS antenna and a transducer on a sea-surface platform.

* `key`: Sea surface plataform numebr to be optimized
* `lat`: Site latitude
* `TR_DEPTH`: Transducer depth from the sea-surface
* `NPB`: Number of temporal B-spline bases
* `eps`: Convergence threshold (`eps=1.e-4` by default)
* `IMAX`: Maximum number of iterations (`IMAX=50` by default)
* `delta_pos`: Infinitesimal amount of the array displacements to calculate the Jacobian matrix (`delta_pos=1.e-4`)
* `delta_offset`: Infinitesimal amount of the offset between a GNSS antenna and a transducer (`delta_offset=1.e-4`)
* `fn1`: Input file name for an offset between a GNSS antenna and a transducer on a sea-surface platform [m] (`fn1="tr-ant.inp"` by default)
* `fn2`: Input file name for the initial seafloor transponder positions [m] (`fn2="pxp-ini.inp"` by default)
* `fn3`: Input file name for the initial sound speed profile (`fn3="ss_prof.inp"` by default)
* `fn4`: Input file name for the basic observational data  (`fn4="obsdata.inp"` by default)
* `fno0`: Output file name for logging  (`fno0=log.txt` by default)
* `fno1`: Output file name for the estimated parameters and their stds (`fno1=solve.out` by default)
* `fno2`: Output file name for the estimated array displacement (`fno2=position.out` by default)
* `fno3`: Output file name for the residuals (`fno3=residual.out` by default)
* `fno4`: Output file name for the estimated B-spline bases (`fno4=bspline.out` by default)
* `fno5`: Output file name for the estimated offset between GNSS antenna and transducer (`fno5=tr-ant.out` by default)

# Example
    static_array_TR(38.0,5.0,100,delta_offset=5.e-4)
"""
function static_array_TR(key, lat, TR_DEPTH::Vector{Float64},NPB=100::Int64; fn1="tr-ant.inp"::String,fn2="pxp-ini.inp"::String,fn3="ss_prof.inp"::String,fn4="obsdata.inp"::String,eps=1.e-4,ITMAX=50::Int64, delta_pos=1.e-4, delta_offset=1.e-4,fno0="log.txt"::String,fno1="solve.out"::String,fno2="position.out"::String,fno3="residual.out"::String,fno4="bspline.out"::String,fno5="tr-ant.out"::String)
  println(stderr," === GNSS-A positioning: static_array_TR  ===")
  nds0 = size(TR_DEPTH)[1]
  # --- Input check
  if NPB < 1
    error(" static_array_TR: NPB must be more than 1")
  end
  # --- Start log
  time1 = now()
  place = pwd()
  open(fno0,"w") do out0 
  println(out0,time1)
  println(out0,"static_array_TR.jl at $place")
  # --- Set parameters
  println(stderr," --- Set parameters")
  NP0 = 6; NC = 18
  dx = delta_pos; dy = delta_pos; dz = delta_pos; de = delta_offset
  println(out0,"Convergence_eps: $eps")
  println(out0,"Number_of_B-spline_knots: $NPB")
  println(out0,"Default_latitude: $lat")
  println(out0,"Maximum_iterations: $ITMAX")
  println(out0,"Delat_position: $delta_pos")
  println(out0,"Delat_offset: $delta_offset")
  for n in 1:nds0
    println(out0,"TR_DEPTH-$n: $TR_DEPTH[$n]")
  end
  TR_DEPTH0 = minimum(TR_DEPTH)
  # --- Read data
  println(stderr," --- Read files")
  e = read_ant(fn1)
  numk, px, py, pz = read_pxppos(fn2)
  z, v, nz_st, numz = read_prof(fn3,TR_DEPTH0)
  num, nk, tp, t1, x1, y1, z1, h1, p1, r1, t2, x2, y2, z2, h2, p2, r2, nf, ids = read_obsdata(fn4)
  if z[end] < maximum(abs.(pz))                                                 
    error(" static_array_TR: maximum water depth of $fn3 must be deeper than site depth of $fn2")
  end
  nds = size(e)[2]
  println(out0,"Number_of_sea-surface-platforms: $nds")
  # --- Fixed Earth radius
  Rg, Rl = localradius(lat)

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
  d = zeros(num); H = zeros(num,NP); a0 = zeros(NP); a = zeros(NP)
  dc = zeros(num); dr = zeros(num); delta = 1.e6; rms = 1.e6
  sigma2 = 0.0; aic = 0.0; bic = 0.0; Hinv = zeros(NP,NP)
  a0[4] = e[1,key]; a0[5] = e[2,key]; a0[6] = e[3,key]

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
      # --- Set TR pos
      e0 = zeros(3)
      e0[1] = a0[4]; e0[2] = a0[5]; e0[3] = a0[6]
      if ids[n] == key
        xd1[n], yd1[n], zd1[n] = anttena2tr(x1[n],y1[n],z1[n],h1[n],p1[n],r1[n],e0)
        xd2[n], yd2[n], zd2[n] = anttena2tr(x2[n],y2[n],z2[n],h2[n],p2[n],r2[n],e0)
      else
        xd1[n], yd1[n], zd1[n] = anttena2tr(x1[n],y1[n],z1[n],h1[n],p1[n],r1[n],e[:,ids[n]])
        xd2[n], yd2[n], zd2[n] = anttena2tr(x2[n],y2[n],z2[n],h2[n],p2[n],r2[n],e[:,ids[n]])
      end
      # --- Calculate TT
      tc1, to1, vert1 = xyz2tt_rapid(px[k]+a0[1],py[k]+a0[2],pz[k]+a0[3],xd1[n],yd1[n],zd1[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
      tc2, to2, vert2 = xyz2tt_rapid(px[k]+a0[1],py[k]+a0[2],pz[k]+a0[3],xd2[n],yd2[n],zd2[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
      vert = (vert1 + vert2) / 2.0
      tc = tc1 + tc2
      d[n] = (tp[n] - tc)*vert
      # --- Differential for pos
      tcx1, to1, vert1 = xyz2tt_rapid(px[k]+a0[1]+dx,py[k]+a0[2],pz[k]+a0[3],xd1[n],yd1[n],zd1[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
      tcx2, to2, vert2 = xyz2tt_rapid(px[k]+a0[1]+dx,py[k]+a0[2],pz[k]+a0[3],xd2[n],yd2[n],zd2[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
      tcx = tcx1 + tcx2
      tcy1, to1, vert1 = xyz2tt_rapid(px[k]+a0[1],py[k]+a0[2]+dy,pz[k]+a0[3],xd1[n],yd1[n],zd1[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
      tcy2, to2, vert2 = xyz2tt_rapid(px[k]+a0[1],py[k]+a0[2]+dy,pz[k]+a0[3],xd2[n],yd2[n],zd2[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
      tcy = tcy1 + tcy2
      tcz1, to1, vert1 = xyz2tt_rapid(px[k]+a0[1],py[k]+a0[2],pz[k]+a0[3]+dz,xd1[n],yd1[n],zd1[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
      tcz2, to2, vert2 = xyz2tt_rapid(px[k]+a0[1],py[k]+a0[2],pz[k]+a0[3]+dz,xd2[n],yd2[n],zd2[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
      tcz = tcz1 + tcz2
      # --- Difference for tr-ant
      e0[1] = a0[4] + de; e0[2] = a0[5]; e0[3] = a0[6]
      if ids[n] == key
        xd1[n], yd1[n], zd1[n] = anttena2tr(x1[n],y1[n],z1[n],h1[n],p1[n],r1[n],e0)
        xd2[n], yd2[n], zd2[n] = anttena2tr(x2[n],y2[n],z2[n],h2[n],p2[n],r2[n],e0)
      else
        xd1[n], yd1[n], zd1[n] = anttena2tr(x1[n],y1[n],z1[n],h1[n],p1[n],r1[n],e[:,ids[n]])
        xd2[n], yd2[n], zd2[n] = anttena2tr(x2[n],y2[n],z2[n],h2[n],p2[n],r2[n],e[:,ids[n]])
      end
      tc1, to1, vert1 = xyz2tt_rapid(px[k]+a0[1],py[k]+a0[2],pz[k]+a0[3],xd1[n],yd1[n],zd1[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
      tc2, to2, vert2 = xyz2tt_rapid(px[k]+a0[1],py[k]+a0[2],pz[k]+a0[3],xd2[n],yd2[n],zd2[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
      tce1 = tc1 + tc2
      e0[1] = a0[4]; e0[2] = a0[5] + de; e0[3] = a0[6]
      if ids[n] == key
        xd1[n], yd1[n], zd1[n] = anttena2tr(x1[n],y1[n],z1[n],h1[n],p1[n],r1[n],e0)
        xd2[n], yd2[n], zd2[n] = anttena2tr(x2[n],y2[n],z2[n],h2[n],p2[n],r2[n],e0)
      else
        xd1[n], yd1[n], zd1[n] = anttena2tr(x1[n],y1[n],z1[n],h1[n],p1[n],r1[n],e[:,ids[n]])
        xd2[n], yd2[n], zd2[n] = anttena2tr(x2[n],y2[n],z2[n],h2[n],p2[n],r2[n],e[:,ids[n]])
      end
      tc1, to1, vert1 = xyz2tt_rapid(px[k]+a0[1],py[k]+a0[2],pz[k]+a0[3],xd1[n],yd1[n],zd1[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
      tc2, to2, vert2 = xyz2tt_rapid(px[k]+a0[1],py[k]+a0[2],pz[k]+a0[3],xd2[n],yd2[n],zd2[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
      tce2 = tc1 + tc2
      e0[1] = a0[4]; e0[2] = a0[5]; e0[3] = a0[6] + de
      if ids[n] == key
        xd1[n], yd1[n], zd1[n] = anttena2tr(x1[n],y1[n],z1[n],h1[n],p1[n],r1[n],e0)
        xd2[n], yd2[n], zd2[n] = anttena2tr(x2[n],y2[n],z2[n],h2[n],p2[n],r2[n],e0)
      else
        xd1[n], yd1[n], zd1[n] = anttena2tr(x1[n],y1[n],z1[n],h1[n],p1[n],r1[n],e[:,ids[n]])
        xd2[n], yd2[n], zd2[n] = anttena2tr(x2[n],y2[n],z2[n],h2[n],p2[n],r2[n],e[:,ids[n]])
      end
      tc1, to1, vert1 = xyz2tt_rapid(px[k]+a0[1],py[k]+a0[2],pz[k]+a0[3],xd1[n],yd1[n],zd1[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
      tc2, to2, vert2 = xyz2tt_rapid(px[k]+a0[1],py[k]+a0[2],pz[k]+a0[3],xd2[n],yd2[n],zd2[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
      tce3 = tc1 + tc2
      # --- Fill matrix
      H[n,1] = (tcx-tc)/dx*vert; H[n,2]=(tcy-tc)/dy*vert; H[n,3]=(tcz-tc)/dz*vert
      H[n,4] = (tce1-tc)/de*vert; H[n,5]=(tce2-tc)/de*vert; H[n,6]=(tce3-tc)/de*vert
      if it == 1
        for m in 1:NPB
          if id[m] >= 1
            b0 = zeros(NPB)
            b0[m] = 1.0
            H[n,NP0+id[m]] = tbspline3((t1[n]+t2[n])/2.0,ds,tb,b0,NPB)
          end
        end
      end
    end
    Hinv = inv(transpose(H)*H)
    a = Hinv*transpose(H)*d
    dc = H*a
    dr = d - dc
    rms = std(dr)
    sa = num * rms^2
    sigma2 = sa / (num - NP)
    delta = std(a[1:NP0])
    a0[1:NP0] += a[1:NP0]
    a0[NP0+1:NP] = a[NP0+1:NP]
    println(stderr," Temporal position: $(a0[1:3]), $delta, $rms")
    println(stderr," Temporal offset: $(a0[4:6])")
    println(out0,"     Iteration: $it $(a0[1]) $(a0[2]) $(a0[3]) $delta $rms")
    println(out0,"                 $(a0[4]) $(a0[5]) $(a0[6])")
    it += 1
  end
  println(stderr," End of loop ",it-1)
  println(out0,"End of iteration")
  cv = sqrt.(sigma2*abs.(diag(Hinv))) # Error
  tc = (t1[1] + t2[num]) / 2.0
  a = transpose(a0[1:3])
  # --- Fill NTD basis
  b = zeros(NPB)
  for m in 1:NPB
    if id[m] >= 1
      b[m] = a0[NP0+id[m]]
    else
      b[m] = 0.0
    end
  end
  td = zeros(num)
  for n in 1:num
    td[n] = tbspline3((t1[n]+t2[n])/2.0,ds,tb,b,NPB)
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
    tt = (t1 + t2) / 2.0
    for i in 1:num
      println(out,"$(tt[i]) $(nk[i]) $(d[i]) $(dc[i]) $(dr[i])")
    end
  end
  open(fno4,"w") do out
    for i in 1:NPB                                                              
      println(out,"$i $(id[i]) $(b[i])")
    end
  end
  open(fno5,"w") do out
    println(out,"# TR_offset: dx, dy, dz, sigmax, sigmay, sigmaz")
    println(out,"$(a0[4]) $(a0[5]) $(a0[6]) $(cv[4]) $(cv[5]) $(cv[6])")
  end

# --- Close process --- #
  time2 = now()
  println(stderr," Start time:",time1)
  println(stderr," Finish time:",time2)
  println(out0,time2)
  end
end
