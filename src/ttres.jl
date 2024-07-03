#using Dates
#using Statistics
#using LinearAlgebra

export ttres
"""
    ttres(lat,TR_DEPTH; fn1,fn2,fn3,fn4,fno,fno0,save)

Calculate travel-time residuals using fixed transponder positions.

Input:
* `lat`: Site latitude
* `TR_DEPTH`: Transducer depth from the sea-surface
* `NPB`: Number of temporal B-spline bases
* `fn1`: Input file name for an offset between a GNSS antenna and a transducer on a sea-surface platform [m] (`fn1="tr-ant.inp"` by default)
* `fn2`: Input file name for the initial seafloor transponder positions [m] (`fn2="pxp-ini.inp"` by default)
* `fn3`: Input file name for the initial sound speed profile (`fn3="ss_prof.inp"` by default)
* `fn4`: Input file name for the basic observational data  (`fn4="obsdata.inp"` by default)
* `fno0`: Output file name for logging  (`fno0=log.txt` by default)
* `fno`: Output file name for the travel-time residuals (`fno=ttres.out` by default)
* `save`: if `save=true`, calculation results are saved in `fno`

Output:
* `nv`: Shot number
* `kv`: Seafloor transponder number
* `t1`: Signal transmitting time [sec]
* `t2`: Signal recieving time [sec]
* `tp`: Observed travel-time
* `tc`: Calculated travel-time
* `tr`: Travel-time residual
* `vert`: Normalizing factor

# Example
    nv, kv, t1, t2, tp, tc, tr, vert = ttres(lat,TR_DEPTH)
"""
function ttres(lat=38.0,TR_DEPTH=[3.0];fn1="tr-ant.inp", fn2="pxp-ini.inp", fn3="ss_prof.inp", fn4="obsdata.inp",fno="ttres.out",fno0="log.txt",save=false)
  println(stderr," === Calculate trave-time residuals  ===")
  # --- Start log
  time1 = now()
  place = pwd()
  open(fno0,"w") do out0 
  println(out0,time1)
  println(out0,"ttres.jl at $place")
  # --- Set parameters
  println(stderr," --- Set parameters")
  nds0 = size(TR_DEPTH)[1]
  NC = 18 # Number of fixed parameters
  println(out0,"Default_latitude: $lat")
  println(out0,"TR_DEPTH: $TR_DEPTH")
  println(out0,"Transducer_Antenna: $fn1")
  println(out0,"Transponder_position: $fn2")
  println(out0,"Sound_speed_profile: $fn3")
  println(out0,"Observational_data: $fn4")
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
    error(" ttres: maximum water depth of $fn3 must be deeper than site depth of $fn2")
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

# --- Main Anlysis --- #
  nv = Int64[]; kv = Int64[]; tcv = Float64[]; vertv = Float64[]
  if save == true
    open(fno,"w") do out
      for n in 1:num
        k = nk[n]  # PXP number
        Rg, Rl = localradius(lat)
        # --- Calculate TT
        tc1, Nint1, vert1 = xyz2tt(px[k],py[k],pz[k],xd1[n],yd1[n],zd1[n],z,v,nz_st,numz,Rg,TR_DEPTH[ids[n]])
        tc2, Nint2, vert2 = xyz2tt(px[k],py[k],pz[k],xd2[n],yd2[n],zd2[n],z,v,nz_st,numz,Rg,TR_DEPTH[ids[n]])
        vert = (vert1 + vert2) / 2.0
        tc = tc1 + tc2
        println(out,"$n $k $(t1[n]) $(t2[n]) $(tp[n]) $tc $(tp[n]-tc) $vert")
        push!(nv,n); push!(kv,k); push!(tcv,tc); push!(vertv,vert)
      end
    end
  else
    for n in 1:num
      k = nk[n]  # PXP number
      Rg, Rl = localradius(lat)
      # --- Calculate TT
      tc1, Nint1, vert1 = xyz2tt(px[k],py[k],pz[k],xd1[n],yd1[n],zd1[n],z,v,nz_st,numz,Rg,TR_DEPTH[ids[n]])
      tc2, Nint2, vert2 = xyz2tt(px[k],py[k],pz[k],xd2[n],yd2[n],zd2[n],z,v,nz_st,numz,Rg,TR_DEPTH[ids[n]])
      vert = (vert1 + vert2) / 2.0
      tc = tc1 + tc2
      push!(nv,n); push!(kv,k); push!(tcv,tc); push!(vertv,vert)
    end
  end
  return nv, kv, t1, t2, tp, tcv, tp-tcv, vertv


# --- Close process --- #
  time2 = now()
  println(stderr," Start time:",time1)
  println(stderr," Finish time:",time2)
  println(out0,time2)
  end
end
