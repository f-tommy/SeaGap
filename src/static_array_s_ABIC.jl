#using Dates
#using Statistics
#using LinearAlgebra

export static_array_s_ABIC
"""
    static_array_s_ABIC(r1,r2,r3,lat,TR_DEPTH,NPB1,NPB2; fn1,fn2,fn3,fn4,eps,ITMAX,delta_pos,fno0,fno1,fno2,fno3,fno4,fno5,fno6)

Perform static array positioning considering shallow gradients using tranducer-based observational data file.
Long-period and short-period NTDs are modeled by small (NPB1) and large (NPB2) numbers of 3d B-spline bases.
Range for the hyper-parameter constraining the norm of S-NTD (alpha) to be investigated is given by (`r1`,`r2`) with interval of `r3`.

* `lat`: Site latitude
* `TR_DEPTH`: Transducer depth from the sea-surface
* `NPB1`: Number of L-NTD B-spline bases
* `NPB2`: Number of S-NTD B-spline bases
* `eps`: Convergence threshold (`eps=1.e-4` by default)
* `IMAX`: Maximum number of iterations (`IMAX=50` by default)
* `delta_pos`: Infinitesimal amount of the array displacements to calculate the Jacobian matrix (`delta_pos=1.e-4`)
* `fn1`: Input file name for the antenna-TR offset (`fn1="tr-ant.inp"` by default)
* `fn2`: Input file name for the initial seafloor transponder positions [m] (`fn2="pxp-ini.inp"` by default)
* `fn3`: Input file name for the initial sound speed profile (`fn3="ss_prof.inp"` by default)
* `fn4`: Input file name for the observational data (`fn1="obsdata.inp"` by default)
* `fno0`: Output file name for logging  (`fno0=log.txt` by default)
* `fno`: Output file name for the ABIC (`fno=ABIC_search.out` by default)

# Example
    static_array_s_ABIC(-2.0,2.0,0.5,36.2,[0.8],5,100)
"""
function static_array_s_ABIC(r1,r2,r3,lat,TR_DEPTH::Vector{Float64},NPB1=5,NPB2=100::Int64; fn1="tr-ant.inp"::String,fn2="pxp-ini.inp"::String,fn3="ss_prof.inp"::String,fn4="obsdata.inp"::String,eps=1.e-4,ITMAX=50::Int64, delta_pos=1.e-4, fno0="log.txt"::String,fno="ABIC_search.out"::String,delete=false::Bool,spc=false)
  println(stderr," === GNSS-A positioning: static_array_s_ABIC  ===")
  # --- Input check
  nds0 = size(TR_DEPTH)[1]
  if NPB1 < 2
    error(" static_array_s_ABIC: NPB1 must be more than 2")
  end
  if NPB2 < 2
    error(" static_array_s_ABIC: NPB2 must be more than 2")
  end
  if r1 > r2
    error(" static_array_s_ABIC: r2 must be larger than r1")
  end
  range = r1:r3:r2
  # --- Start log
  time1 = now()
  place = pwd()
  open(fno0,"w") do out0 
  println(out0,time1)
  println(out0,"static_array_s_ABIC.jl at $place")
  # --- Set parameters
  println(stderr," --- Set parameters")
  NP0 = 5; NC = 18 # Number of fixed parameters
  dx = delta_pos; dy = delta_pos; dz = delta_pos
  println(out0,"Convergence_eps: $eps")
  println(out0,"Hyper-parameter_for_S-NTD: 10^ $alpha")
  println(out0,"Number_of_B-spline_knots1: $NPB1")
  println(out0,"Number_of_B-spline_knots2: $NPB2")
  println(out0,"Default_latitude: $lat")
  println(out0,"Maximum_iterations: $ITMAX")
  println(out0,"Delta_position: $delta_pos")
  println(out0,"TR_DEPTH: $TR_DEPTH")
  for n in 1:nds0
    println(out0,"TR_DEPTH-$n: $TR_DEPTH[$n]")
  end
  TR_DEPTH0 = minimum(TR_DEPTH)
  # --- Delete pre-existence file "ABIC_search.out"
  println(out0,"Output file: $fno")
  if delete == true
    rm(fno)
    println(out0,"Remove the old version of $fno")
  end
  # --- Read data
  println(stderr," --- Read files")
  e = read_ant(fn1)
  numk, px, py, pz = read_pxppos(fn2)
  z, v, nz_st, numz = read_prof(fn3,TR_DEPTH0)
  num, nk, tp, t1, x1, y1, z1, h1, p1, r1, t2, x2, y2, z2, h2, p2, r2, nf, ids = read_obsdata(fn4)
  if z[end] < maximum(abs.(pz))                                                 
    error(" static_array: maximum water depth of $fn3 must be deeper than site depth of $fn2")
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
  println(stderr,"     Transducer_height:",tr_height)
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
  NPBV2, id2 = retrieveb(NPB2,tb2,ds2,t1,t2,num) 
  NPA = 1
  if spc == false
    NPBV2 = zeros(Int64,1)
    id2 = zeros(Int64,NPB2,1)
    NPBV2[1], id2[:,1] = retrieveb(NPB2,tb2,ds2,t1,t2,num) 
  else
    NPBV2 = zeros(Int64,ntr)
    id2 = zeros(Int64,NPB2,ntr)
    NPA = ntr
    for k in 1:ntr
      t1e = t1[ids .== k]
      t2e = t2[ids .== k]
      nume = size(t1e)[1]
      NPBV2[k], id2[:,k] = retrieveb(NPB2,tb2,ds2,t1e,t2e,nume) 
    end
  end 
  # --- Initialize
  NP = NP0 + NPBV1 + sum(NPBV2)
  d = zeros(num); H = zeros(num,NP); a0 = zeros(NP); a = zeros(NP)
  dc = zeros(num); dr = zeros(num); delta = 1.e6; rms = 1.e6
  G0 = zeros(NP,NP)
  d1 = zeros(num); d2 = zeros(num); d3 = zeros(num)
  xd = zeros(num); yd = zeros(num)
  sigma2 = 0.0; abic = 0.0; Hinv = zeros(NP,NP)
  Rg, Rl = localradius(lat)

# --- Main Anlysis --- #
  for alpha in range
    open(fno,"a") do outa
    # --- Initialize
    NP = NP0 + NPBV1 + sum(NPBV2)
    d = zeros(num); H = zeros(num,NP); a0 = zeros(NP); a = zeros(NP)
    dc = zeros(num); dr = zeros(num); delta = 1.e6; rms = 1.e6
    G0 = zeros(NP,NP)
    d1 = zeros(num); d2 = zeros(num); d3 = zeros(num)
    xd = zeros(num); yd = zeros(num)
    sigma2 = 0.0; abic = 0.0; Hinv = zeros(NP,NP)
    for i in NP0+NPBV1+1:NP
      G0[i,i] = 10.0^alpha
    end

    # --- Inversion
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
        tc1, to1, vert1 = xyz2tt_rapid(px[k]+a0[1],py[k]+a0[2],pz[k]+a0[3],xd1[n],yd1[n],zd1[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
        tc2, to2, vert2 = xyz2tt_rapid(px[k]+a0[1],py[k]+a0[2],pz[k]+a0[3],xd2[n],yd2[n],zd2[n],Rg,Tv0[k],Vd[k],Vr[k],tr_height,cc[k,1:NC])
        vert = (vert1 + vert2) / 2.0
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
        H[n,4] = xd[n]; H[n,5] = yd[n]
        if it == 1
          for m in 1:NPB1
            if id1[m] >= 1
              b0 = zeros(NPB1)
              b0[m] = 1.0
              H[n,NP0+id1[m]] = tbspline3((t1[n]+t2[n])/2.0,ds1,tb1,b0,NPB1)
            end
          end
          if spc == false
            ll = 1
          else
            ll = ids[n]
          end
          for m in 1:NPB2
            if id2[m,ll] >= 1
              b0 = zeros(NPB2)
              b0[m] = 1.0
              if ll == 1
                H[n,NP0+NPBV1+id2[m,ll]] = tbspline3((t1[n]+t2[n])/2.0,ds2,tb2,b0,NPB2)
              else
                H[n,NP0+NPBV1+sum(NPBV2[1:ll-1])+id2[m,ll]] = tbspline3((t1[n]+t2[n])/2.0,ds2,tb2,b0,NPB2)
              end
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
      d1 = H*a1; d2 = H*a2; d3 = H*a3
      rms = std(dr)
      sa = LinearAlgebra.dot(dr,dr)
      sigma2 = sa / (num-NP)
      abic = num*log(sa) - sum(NPBV2)*log(10.0^alpha) + logdet(Horg)
      delta = std(a[1:3])
      a0[1:3] += a[1:3]
      a0[4:NP] = a[4:NP]
      println(stderr," Temporal position: $(a0[1:3]), $delta, $rms")
      println(out0,"     Iteration: $it $(a0[1]) $(a0[2]) $(a0[3]) $delta $rms")
      it += 1
    end
    println(stderr," End of loop ",it-1)
    println(out0," End of loop for $alpha",it-1)
    println(outa,"$alpha $NPB1 $NPB2 0.0 $abic")
  end
  end

# --- Close process --- #
  time2 = now()
  println(stderr," Start time:",time1)
  println(stderr," Finish time:",time2)
  println(out0,time2)
  end
end
