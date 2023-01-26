#using DelimitedFiles
#using Statistics
#using Dates
#using Printf
#include("unixsort.jl")
#include("dateprocessing.jl")
# Usage: convert_displacement(1.79,-1.79)

export bitarr_to_int, convert_displacement

function bitarr_to_int(arr)
  return sum(arr .* (2 .^ collect(length(arr)-1:-1:0)))
end

"""
    convert_displacement(vx,vy,vz; redu,sredu_hor,sredu_ver,fno,fn,t0)

Convert a position time-series file `fn` (1: the observation period [sec], 2-4: EW-NS-UD positions [m],5-7: EW-NS-UD Std of positions [m]) into an easily-readable time-series file `fno` (1: the observation period [year], 2-4: EW-NS-UD position [cm] eliminating a stable motion `vx`-`vy`-`vz`, 5-7: EW-NS-UD Std of positions [cm], 8-10: ID for usable/unusable data (1/0) for each component).

* `vx`, `vy`, `vz`: constant velocities which are eliminated from the time-series (Unit: cm/yr)
* `redu`: if `redu=true`, the large estimation errors are identified as unusable (estimation errors in EW and NS components > `sredu_hor`; those in UD component > `sredu_ver`); if `redu=false`, all are identified as usable
* `sredu_hor`: Estimation error limit to be usable for horizontal components (`sredu_hor=30.0` [cm] in default)
* `sredu_ver`: Estimation error limit to be usable for vertical component (`sredu_ver=30.0` [cm] in default)
* `fno`: Arranged time-series file (`fno="converted_position.out"` in defalult)
* `fn`: Input time-series file (`fn="position_merge.out"` in defalult)
* `t0`: Reference time (`t0="2000-01-01T12:00:00"` in default; refer Date & Time processing [e.g., `sec2year()`])

# Example
    convert_displacement(1.79,-1.79,0.0,sredu_hor=30.0,sredu_ver=30.0,fno="converted_position.out")
"""
function convert_displacement(vx0=0.0,vy0=0.0,vz0=0.0; redu=true,sredu_hor=30.0,sredu_ver=30.0,fno="converted_position.out"::String,fn="position_merge.out"::String,t0="2000-01-01T12:00:00")
  # --- Read
  println(stderr," --- Read $fn")
  dat0 = DelimitedFiles.readdlm(fn)
  dat = unixsort(dat0,1)
  num, mm = size(dat)
  if num < 1
    error(" No data is given in $fn")
  end
  if mm < 7
    error(" Insufficient data: Need (j2000time DispX DispY DispZ SigmaX SigmaY SigmaZ)")
  end
  j = dat[:,1]; dx0 = dat[:,2]; dy0 = dat[:,3]; dz0 = dat[:,4]; sx0 = dat[:,5]; sy0 = dat[:,6]; sz0 = dat[:,7]
  t = sec2year.(j,t0)
  # --- Remove stable motions
  println(stderr," --- Remove stable motion: x $vx0, y $vy0, z $vz0")
  dx = dx0*100 - vx0 * (t .- t[1])
  dy = dy0*100 - vy0 * (t .- t[1])
  dz = dz0*100 - vz0 * (t .- t[1])
  sx = sx0*100; sy = sy0*100; sz = sz0*100
  dat0 = hcat(t,dx,dy,dz,sx,sy,sz)
  # --- Identify outlier
  println(stderr," --- Identify outliers: $sredu_hor, $sredu_ver")
  if redu == true
    idx = bitarr_to_int([dat0[:,5] .< sredu_hor])
    idy = bitarr_to_int([dat0[:,6] .< sredu_hor])
    idz = bitarr_to_int([dat0[:,7] .< sredu_ver])
  else
    idx = ones(Int64,num)
    idy = ones(Int64,num)
    idz = ones(Int64,num)
  end
  open(fno,"w") do out
    println(out,"# Year DispX DispY DispZ SigmaX SigmaY SigmaZ UsableX UsableY UsableZ (detrend $vx0 $vy0 $vz0)")
    for n in 1:num
      @printf(out,"%4.5f %3.3f %3.3f %3.3f %3.3f %3.3f %3.3f %d %d %d\n",dat0[n,1],dat0[n,2],dat0[n,3],dat0[n,4],dat0[n,5],dat0[n,6],dat0[n,7],idx[n],idy[n],idz[n])
    end
  end
  println(stderr," Output converted position: $fno")
end
