#using Dates
#using Printf

export obsdata_format
"""
    obsdata_format(numk; fno0,fno,fn1,fn21,fn22,maxp)

Make observational data file `fno` from GNSS antenna position time-series with attitude `fn1` and travel-time time-series file `fn21`-k.`fn22`.

* `numk`: Number of seafloor tranponders at a site
* `fn0`: Log file name
* `fno`: Output file name (`fno="obsdata.inp"` in default)
* `fn1`: File name for GNSS antenna position time-series with attitude (`fn1=gps.jxyhhpr` in default)
* `fn21` and `fn22`: File names for travel-time time-series file (`fn21="pxp-"` and `fn22=".jttq"` in default; if default names with 3 tranponders, "pxp-1.jttq", "pxp-2.jttq", and "pxp-3.jttq" must be prepared)
* `maxp`: Maximum number for readable shots (`maxp=50000` in default; if you have shot data > 50000, you have to provide larger value for `maxp`)

# Example
    obsdata_format(4)
"""
function obsdata_format(numk=4::Int64 ;fno0="log.txt"::String,fno="obsdata.inp"::String,fn1="gps.jxyhhpr"::String,fn21="pxp-"::String,fn22=".jttq"::String,maxp = 50000::Int64)
  println(stderr," === GNSS-A observational data formatting  ===")
  # --- Log
  time1 = now()
  place = pwd()
  open(fno0,"w") do out0
  println(out0,time1)
  println(out0,"obsdata_format.jl at $place")
  # --- Set parameters
  println(stderr," --- Set parameters")
  # --- Read data
  println(stderr," --- Read files")
  numj, tg0, xg0, yg0, zg0, hg0, pg0, rg0 = read_gps(fn1)
  tp1 = zeros(numk,maxp); tt0 = zeros(numk,maxp); tt = zeros(numk,maxp); qq = zeros(numk,maxp)
  nump = zeros(Int64,numk)
  for k in 1:numk
    fn = fn21*"$k"*fn22
    nump[k] = read_jttq(fn,k,tp1,tt0,tt,qq)
  end

# --- Formatting --- #
  # --- Interpolate gps position
  println(stderr," --- Interpolate when transmitting and recieving")
  xg1,yg1,zg1,hg1,pg1,rg1,xg2,yg2,zg2,hg2,pg2,rg2,tp2 = interpolate_gps(tg0,xg0,yg0,zg0,hg0,pg0,rg0,tp1,tt0,nump,numk,maxp)
  # --- Check simultaneous shot
  println(stderr," --- Check simultaneous shot")
  t_all = []
  for k in 1:numk
    append!(t_all,tp1[k,1:nump[k]])
  end
  tp = sort(unique(t_all))
  unum = size(tp)[1]

# --- Output --- #
  dat0 = []
  for k in 1:numk
    kk = nump[k]
    for m in 1:kk
      ip = 0
      for i in 1:unum
        if tp[i] == tp1[k,m]
          ip = i
        end
      end
      push!(dat0,vcat(k,tt[k,m],tp1[k,m],xg1[k,m],yg1[k,m],zg1[k,m],hg1[k,m],pg1[k,m],rg1[k,m],tp2[k,m],xg2[k,m],yg2[k,m],zg2[k,m],hg2[k,m],pg2[k,m],rg2[k,m],ip))
    end
  end
  dat = transpose(hcat(dat0...))
  a = unixsort2(dat,3,1)
  num = size(dat)[1]
  
  open(fno,"w") do out3
  for n in 1:num
    @printf(out3,"%i %2.6f %10.6f %6.6f %6.6f %6.6f %5.5f %3.5f %3.5f %10.6f %6.6f %6.6f %6.6f %5.5f %3.5f %3.5f %d\n",a[n,1],a[n,2],a[n,3],a[n,4],a[n,5],a[n,6],a[n,7],a[n,8],a[n,9],a[n,10],a[n,11],a[n,12],a[n,13],a[n,14],a[n,15],a[n,16],a[n,17])
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
