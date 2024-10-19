#using Dates
#using Statistics
#using LinearAlgebra

export static_array_AICBIC
"""
    static_array_AICBIC(r1,r2,r3,lat,TR_DEPTH; fn1,fn2,fn3,fn4,eps,ITMAX,delta_pos,delete,fno0,fno)

Perform static array positioning and calculate AIC and BIC values for various number of temporal B-spline bases.
Range for the number of temporal B-spline bases to be investigated is given by (`r1`,`r2`) with interval of `r3`.

* `lat`: Site latitude
* `TR_DEPTH`: Transducer depth from the sea-surface
* `eps`: Convergence threshold (`eps=1.e-4` by default)
* `IMAX`: Maximum number of iterations (`IMAX=50` by default)
* `delta_pos`: Infinitesimal amount of the array displacements to calculate the Jacobian matrix (`delta_pos=1.e-4`)
* `fn1`: Input file name for an offset between a GNSS antenna and a transducer on a sea-surface platform [m] (`fn1="tr-ant.inp"` by default)
* `fn2`: Input file name for the initial seafloor transponder positions [m] (`fn2="pxp-ini.inp"` by default)
* `fn3`: Input file name for the initial sound speed profile (`fn3="ss_prof.inp"` by default)
* `fn4`: Input file name for the basic observational data  (`fn4="obsdata.inp"` by default)
* `fno0`: Output file name for logging (`fno0=log.txt` by default)
* `fno`: Output file name for AIC and BIC values (`fno=AICBIC_search.out`)

# Example
    static_array_AICBIC(30,150,5,38.1,2.0)
"""
function static_array_AICBIC(r1::Int64,r2::Int64,r3::Int64,lat,TR_DEPTH::Vector{Float64}; eps=1.e-4,ITMAX=50::Int64, delta_pos=1.e-4, delete=false::Bool, fn1="tr-ant.inp"::String, fn2="pxp-ini.inp"::String, fn3="ss_prof.inp"::String, fn4="obsdata.inp"::String, fno="AICBIC_search.out"::String,fno0="log.txt"::String,spc=false)
  println(stderr," === GNSS-A positioning: static_array_AICBIC  ===")
  # --- Input check
  nds0 = size(TR_DEPTH)[1]
  if r1 > r2
    error(" static_array_AICBIC: r2 must be larger than r1")
  end
  range = r1:r3:r2
  # --- Start log
  time1 = now()
  place = pwd()
  open(fno0,"w") do out0 
  println(out0,time1)
  println(out0,"static_array_AICBIC.jl at $place")
  # --- Delete pre-existence file "AICBIC_search.out"
  println(out0,"Output file: $fno")
  if delete == true
    rm(fno)
    println(out0,"Remove the old version of $fno")
  end
  # --- Set parameters
  println(stderr," --- Set parameters")
  NP0 = 3; NC = 18; dx = delta_pos; dy = delta_pos; dz = delta_pos
  ntr = size(TR_DEPTH)[1]
  println(out0,"Convergence eps: $eps")
  println(out0,"Default latitude: $lat")
  println(out0,"Maximum_iterations: $ITMAX")
  println(out0,"Delta_position: $delta_pos")
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
    error(" static_array_AICBIC: maximum water depth of $fn3 must be deeper than site depth of $fn2")
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
  Rg, Rl = localradius(lat)

# --- Loop for optimal NTD basis knot search --- #
  for NPB in range
    open(fno,"a") do outa
    # --- Set B-spline function
    println(stderr," --- NTD basis $NPB")
    smin, smax, ds, tb = mktbasis(NPB,t1,t2,num)
    NPA = 1
    if spc == false
      NPBV = zeros(Int64,1)
      id = zeros(Int64,NPB,1)
      NPBV[1], id[:,1] = retrieveb(NPB,tb,ds,t1,t2,num) 
    else
      NPBV = zeros(Int64,ntr)
      id = zeros(Int64,NPB,ntr)
      NPA = ntr
      for k in 1:ntr
        t1e = t1[ids .== k]
        t2e = t2[ids .== k]
        nume = size(t1e)[1]
        NPBV[k], id[:,k] = retrieveb(NPB,tb,ds,t1e,t2e,nume) 
      end
    end
    # --- Initialize
    NP = NP0 + sum(NPBV)
    d = zeros(num); H = zeros(num,NP); a0 = zeros(NP); a = zeros(NP)
    dc = zeros(num); dr = zeros(num); delta = 1.e6; rms = 1.e6
    sigma2 = 0.0; aic = 0.0; bic = 0.0; Hinv = zeros(NP,NP)
    aic = 0.0; bic = 0.0

    # --- Inversion
    it = 1
    while delta > eps
      if it > ITMAX
        break
      end
      println(stderr,"Inversion for $NPB $it")
      println(out0,"Inversion for $NPB $it")
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
        H[n,1] = (tcx-tc)/dx*vert; H[n,2]=(tcy-tc)/dy*vert; H[n,3]=(tcz-tc)/dz*vert
        if it == 1
          if spc == false
            ll = 1
          else
            ll = ids[n]
          end
          for m in 1:NPB
            if id[m,ll] >= 1
              b0 = zeros(NPB)
              b0[m] = 1.0
              if ll == 1
                H[n,NP0+id[m,ll]] = tbspline3((t1[n]+t2[n])/2.0,ds,tb,b0,NPB)
              else
                H[n,NP0+sum(NPBV[1:ll-1])+id[m,ll]] = tbspline3((t1[n]+t2[n])/2.0,ds,tb,b0,NPB)
              end
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
      aic = num*log(sa/num) + NP*2.0
      bic = num*log(sa/num) + NP*log(num)
      delta = std(a[1:NP0])
      a0[1:NP0] += a[1:NP0]
      a0[NP0+1:NP] = a[NP0+1:NP]
      println(stderr," Temporal position: $(a0[1:3]), $delta, $rms")
      println(out0,"     Iteration: $it $(a0[1]) $(a0[2]) $(a0[3]) $delta $rms")
      it += 1
    end
    println(stderr," End of loop ",it-1)
    println(out0,"End inversion for $NPB $(it-1)")
    println(outa,"$NPB $aic $bic $rms")
  end
  end

# --- Close process --- #
  time2 = now()
  println(stderr," Start time:",time1)
  println(stderr," Finish time:",time2)
  println(out0,time2)
  end
end
