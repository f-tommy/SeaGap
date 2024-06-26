#using Dates
#using Statistics
#using LinearAlgebra

export obsant2obstr
"""
    obsant2obstr(;fn1,fn2,fno)

Tranform GNSS antenna position in "obsdata.inp"(`fn`) into transducer position in `fno`

* `fn1`: Input file name for an offset between a GNSS antenna and a transducer on a sea-surface platform [m] (`fn1="tr-ant.inp"` by default)
* `fn2`: Input file name for the basic observational data  (`fn4="obsdata.inp"` by default)
* `fno`: Output file name for new obsdata.inp  (`fno=obsdata_tr.inp` by default)

# Example
    obsant2obstr(fn2="obsdata.inp",fno="obsdata_tr1.inp")
"""
function obsant2obstr(ids=1::Int64; fn1="tr-ant.inp"::String,fn2="obsdata.inp"::String,fno="obsdata_tr.inp"::String)
  println(stderr," === Convert Ant2TR in obsdata: obsant2obstr  ===")
  # --- Read data
  println(stderr," --- Read files")
  e = read_ant(fn1)
  num, nk, tp, t1, x1, y1, z1, h1, p1, r1, t2, x2, y2, z2, h2, p2, r2, nf = read_obsdata(fn2)
  # --- Formatting
  println(stderr," --- Calculate TR positions")
  xd1 = zeros(num); xd2 = zeros(num)
  yd1 = zeros(num); yd2 = zeros(num)
  zd1 = zeros(num); zd2 = zeros(num)
  for i in 1:num
    xd1[i], yd1[i], zd1[i] = anttena2tr(x1[i],y1[i],z1[i],h1[i],p1[i],r1[i],e)
    xd2[i], yd2[i], zd2[i] = anttena2tr(x2[i],y2[i],z2[i],h2[i],p2[i],r2[i],e)
  end
  # --- Output
  open(fno,"w") do out
  for n in 1:num
    @printf(out,"%i %2.6f %10.6f %6.6f %6.6f %6.6f %10.6f %6.6f %6.6f %6.6f %d %d\n",nk[n],tp[n],t1[n],xd1[n],yd1[n],zd1[n],t2[n],xd2[n],yd2[n],zd2[n],nf[n],ids)
  end
  end
end
