#using Dates
#using Statistics
#using LinearAlgebra

export obstr_reformat
"""
    obstr_reformat(dt;fn,fn0)

Shot group numbers are re-assigned with time-interval of `dt` for "obsdata_tr.inp".

* `dt`: Time-interval [sec]
* `fn`: Input file name for the basic observational data  (`fn="obsdata_tr.inp"` by default)
* `fn0`: Copy of `fn`  (`fn0="obsdata_tr.inp_org"` by default)

# Example
    obstr_format(60,fn="obsdata_tr.inp",fn0="obsdata_tr.inp_org")
"""
function obstr_reformat(dt=60; fn="obsdata_tr.inp"::String,fn0="obsdata_tr.inp_org"::String)
  println(stderr," === Reformat shot group number in obsdata_tr.inp  ===")
  # --- Read data
  println(stderr," --- Read files")
  #num, nk, tp, t1, xd1, yd1, zd1, t2, xd2, yd2, zd2, nf, ids = read_obsdata_tr(fn)
  num,np,dat0 = read_matrix(fn)
  dat = unixsort(dat0,3)
  ts0 = floor(dat[1,3]) - 1
  te0 = dat[end,3] + dt
  # --- Formatting
  cp("$fn","$fn0",force=true)
  open(fn,"w") do out
    ts = ts0
    nf = 0
    while ts <= te0
      dat0 = dat[dat[:,3].>=ts,:]
      ts += dt
      nf += 1
      datv = dat0[dat0[:,3].<ts,:]
      NN,MM = size(datv)
      for n in 1:NN
        @printf(out,"%i %2.6f %10.6f %6.6f %6.6f %6.6f %10.6f %6.6f %6.6f %6.6f %d %d\n",datv[n,1],datv[n,2],datv[n,3],datv[n,4],datv[n,5],datv[n,6],datv[n,7],datv[n,8],datv[n,9],datv[n,10],nf,datv[n,12])
      end
    end
  end
end
