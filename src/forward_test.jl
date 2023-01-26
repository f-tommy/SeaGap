#using Random
#using Dates
#using Statistics
#using LinearAlgebra
#using DelimitedFiles
#include("read_gnssa.jl")
#include("traveltime.jl")
# Usage: forward_test(38, 3.0, [1000,1000,-3000],txtout=false)

export forward_test
"""
    forward_test(lat,XDUCER_DEPTH,pos; txtout,fn,fno)                                                                 
Forward calculation for trave-times (See Tutorials in the online manual).

* `lat`: Site latitude
* `XDUCER_DEPTH`: Transducer depth from the sea-surface
* `pos`: A transponder position (3 components vector)
* `txtout`: if `txtout=true`, calculation results are written in `fno`
* `fn`: Input sound speed structure file
* `fno`: Output text file

# Example
    forward_test(38, 3.0, [1000,1000,-3000],txtout=false)
"""
function forward_test(lat, XDUCER_DEPTH, pos::Vector{}; txtout=false::Bool,fn="ss_prof.zv"::String,fno="synthetic.txt"::String)
# --- Set transponder position
  px = pos[1]; py = pos[2]; pz = pos[3]
  z, v, nz_st, numz = read_prof(fn,XDUCER_DEPTH)
  Rg, Rl = localradius(lat)
  if z[end] < abs(pz)                                                 
    error(" forward_test: maximum water depth of $fn3 must be deeper than site depth of $fn2")
  end

# --- Main Anlysis --- #
  xdv = Float64[]; ydv = Float64[]; zdv = Float64[]; tcv = Float64[]
  yd = 0.
  for xd in 0:100:2000
    # --- Set z-value randomly
    zd = rand()*-5.0
    # --- Calculate TT
    tc, Nint, vert = xyz2tt(px,py,pz,xd,yd,zd,z,v,nz_st,numz,Rg,XDUCER_DEPTH)
    push!(xdv,xd); push!(ydv,yd); push!(zdv,zd); push!(tcv,tc)
  end
  if txtout == true
    open(fno,"w") do out
      println(out,"# X, Y, Z, TT: One-way travel-time calculation for transponder($px,$py,$pz) with the profile $fn")
      Base.print_array(out,hcat(xdv,ydv,zdv,tcv))
    end
  end
  return xdv, ydv, zdv, tcv
end
