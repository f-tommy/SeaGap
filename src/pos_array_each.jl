#using Dates
#using Printf
#using Statistics
#using LinearAlgebra

export pos_array_each
"""
    pos_array_each(lat,XDUCER_DEPTH; fn1,fn2,fn3,fn4,NR,eps,ITMAX,delta_pos,fno0,fno1,fno2)

Perform static array positioning with optimizing offset between a GNSS antenna and a transducer on a sea-surface platform.

* `lat`: Site latitude
* `XDUCER_DEPTH`: Transducer depth from the sea-surface
* `NR`: Minimum number of data to perform estimation for each group
* `eps`: Convergence threshold (`eps=1.e-4` by default)
* `ITMAX`: Maximum number of iterations (`ITMAX=20` by default)
* `delta_pos`: Infinitesimal amount of the array displacements to calculate the Jacobian matrix (`delta_pos=1.e-4`)
* `fn1`: Input file name for an offset between a GNSS antenna and a transducer on a sea-surface platform [m] (`fn1="tr-ant.inp"` by default)
* `fn2`: Input file name for the initial seafloor transponder positions [m] (`fn2="pxp-ini.xyh"` by default)
* `fn3`: Input file name for the initial sound speed profile (`fn3="ss_prof.zv"` by default)
* `fn4`: Input file name for the basic observational data  (`fn4="obsdata.inp"` by default)
* `fno0`: Output file name for logging  (`fno0=log.txt` by default)
* `fno1`: Output file name for the ipositioning results (`fno1=array_each.out` by default)
* `fno2`: Output file name for the travel-time residuals (`fno2=residual_each.out` by default)

# Example
    pos_array_each(42.0,4.0)
"""
function pos_array_each(lat,XDUCER_DEPTH=3.0; NR=3::Int64,eps=1.e-4,ITMAX=20::Int64,delta_pos=1.e-4,fn1="tr-ant.inp"::String,fn2="pxp-ini.xyh"::String,fn3="ss_prof.zv"::String,fn4="obsdata.inp"::String,fno1="array_each.out"::String,fno2="residual_each.out"::String,fno0="log.txt"::String)
  println(stderr," === GNSS-A positioning: pos_array_each  ===")
  # --- Input check
  if XDUCER_DEPTH < 0
    error(" pos_array_each: XDUCER_DEPTH must be positive")
  end
  if NR < 3
    error(" pos_array_each: NR must be more than 3")
  end
  # --- Log
  time1 = now()
  place = pwd()
  open(fno0,"w") do out0
  println(out0,time1)
  println(out0,"pos_array_each.jl at $place")
  println(out0,"XDUCER_Depth: $XDUCER_DEPTH, Lat: $lat")
  # --- Set parameters
  println(stderr," --- Set parameters")
  NP = 3       # Number of parameters for each 
  dx = delta_pos; dy = delta_pos; NC = 18; maxp = 50000 # Fixed parameters
  println(out0,"Number of usable responces for each shot: $NR")
  println(out0,"Delta_position: $delta_pos")
  println(out0,"Maximum_number_of_iterations: $ITMAX")
  println(out0,"Convergence rms: $eps")
  Rg, Rl = localradius(lat)
  # --- Read data
  println(stderr," --- Read files")
  e = read_ant(fn1)
  numk, px, py, pz = read_pxppos(fn2)
  z, v, nz_st, numz = read_prof(fn3,XDUCER_DEPTH)
  num, nk, tp, t1, x1, y1, z1, h1, p1, r1, t2, x2, y2, z2, h2, p2, r2, nf = read_obsdata(fn4)
  numf = findmax(nf)[2]
  if z[end] < maximum(abs.(pz))                                                 
    error(" pos_array_each: maximum water depth of $fn3 must be deeper than site depth of $fn2")
  end

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
  println(out0,"Travel-time correction: $NC")
  xducer_height = ( mean(zd1) + mean(zd2) ) / 2.0
  println(stderr,"     xducer_height:",xducer_height)
  Tv0 = zeros(numk); Vd = zeros(numk); Vr = zeros(numk); cc = zeros(numk,NC)
  for k in 1:numk
    Tv0[k], Vd[k], Vr[k], cc[k,1:NC], rms = ttcorrection(px[k],py[k],pz[k],xducer_height,z,v,nz_st,numz,XDUCER_DEPTH,lat)
    println(stderr,"     RMS for PxP-$k: ",rms)
    println(out0,"     RMS for PxP-$k: ",rms)
  end
  
# --- Main analysis --- #
  # --- Open output file
  open(fno1,"w") do out
  open(fno2,"w") do out2
  # --- Start each positioning
  id = zeros(Int64,numf)
  for n in 1:numf
    ip = zeros(Int64,numk)
    i0 = 1
    for i in 1:num
      if nf[i] == n
        id[n] += 1
        ip[i0] = i
        i0 += 1
      end
    end
    if id[n] >= NR
      println(out0,"Positioning for $n $(t1[n])")
      it = 1
      println(stderr," --- Shot $n with $(id[n]) responses")
      # --- Initialize
      a = zeros(NP); a0 = zeros(NP); delta = 1.e3; rms = 1.e3; sigma2 = 1.e3; NN = id[n]
      d = zeros(NN); dc = zeros(NN); dr = zeros(NN)
      H = zeros(NN,NP); Hinv = zeros(NP,NP); cv = zeros(NP)
      ic = zeros(Int64,NN)
      # --- Loop
      while delta > eps
        if it > ITMAX
          break
        end
        if delta > 1.e4
          break
        end
        println(stderr," --- Iteration: $it")
        # --- Set H-matrix
        for i in 1:NN
          ii = ip[i]
          k = nk[ii]
          ic[i] = k
          # --- Calculate TT
          tc1, to1, vert1 = xyz2tt_rapid(px[k]+a0[1],py[k]+a0[2],pz[k],xd1[ii],yd1[ii],zd1[ii],Rg,Tv0[k],Vd[k],Vr[k],xducer_height,cc[k,1:NC])
          tc2, to2, vert2 = xyz2tt_rapid(px[k]+a0[1],py[k]+a0[2],pz[k],xd2[ii],yd2[ii],zd2[ii],Rg,Tv0[k],Vd[k],Vr[k],xducer_height,cc[k,1:NC])
          vert = (vert1 + vert2) / 2.0
          tc = tc1 + tc2
          d[i] = (tp[ii] - tc)*vert
          # --- Differential TT
          tcx1, to1, vert1 = xyz2tt_rapid(px[k]+a0[1]+dx,py[k]+a0[2],pz[k],xd1[ii],yd1[ii],zd1[ii],Rg,Tv0[k],Vd[k],Vr[k],xducer_height,cc[k,1:NC])
          tcx2, to2, vert2 = xyz2tt_rapid(px[k]+a0[1]+dx,py[k]+a0[2],pz[k],xd2[ii],yd2[ii],zd2[ii],Rg,Tv0[k],Vd[k],Vr[k],xducer_height,cc[k,1:NC])
          tcx = tcx1 + tcx2
          tcy1, to1, vert1 = xyz2tt_rapid(px[k]+a0[1],py[k]+a0[2]+dy,pz[k],xd1[ii],yd1[ii],zd1[ii],Rg,Tv0[k],Vd[k],Vr[k],xducer_height,cc[k,1:NC])
          tcy2, to2, vert2 = xyz2tt_rapid(px[k]+a0[1],py[k]+a0[2]+dy,pz[k],xd2[ii],yd2[ii],zd2[ii],Rg,Tv0[k],Vd[k],Vr[k],xducer_height,cc[k,1:NC])
          tcy = tcy1 + tcy2
          # --- Fill matrix
          H[i,1] = (tcx-tc)/dx*vert; H[i,2] = (tcy-tc)/dy*vert; H[i,3] = 1.0
        end
        # --- Inversion
        Hinv = inv(transpose(H)*H)
        a = Hinv*transpose(H)*d
        dc = H*a
        dr = d - dc
        rms = std(dr)
        sa = NN * rms^2
        sigma2 = sa / (NN-NP)
        a0[1:2] += a[1:2]
        a0[3] = a[3]
        delta = std(a[1:2])
        println(stderr," Temporal position: $(a0[1:2]), $delta, $rms")
        println(out0,"    Iteration: $it $(a0[1]) $(a0[2]) $delta $rms")
        it += 1
      end
      println(stderr," End of loops",it-1)
      println(stderr," Hinv:",diag(Hinv))
      println(stderr," --- Final position: $(a0[1:3]), $delta, $rms")
      
      cv = sqrt.(sigma2*abs.(diag(Hinv))) # Error
      println(out,"$(t1[ip[1]]) $(id[n]) $(a0[1]) $(a0[2]) $(a0[3]) $(cv[1]) $(cv[2]) $(cv[3]) $n")
      for i in 1:NN
        println(out2,"$(t1[ip[1]]) $(ic[i]) $(dr[i]) $(it-1) $n $(id[n])")
      end
    end
  end
  end
  end

# --- Close process --- #
  # --- Computational Time
  time2 = now()
  println(out0,time2)
  println(stderr," Start time:",time1)
  println(stderr," Finish time:",time2)
  end
end
