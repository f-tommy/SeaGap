#using Dates
#using Statistics
#using DelimitedFiles
#include("read_gnssa.jl")
#include("simple_inversion.jl")
# Usage :make_initial()

export make_initial
"""
    make_initial(;fn1,fn2,fn3,fno,error_scale,tscale)

Make an initial data file `fno` for `static_array_mcmcgrad()` or `static_array_mcmcgradc()` using initial transponder positions `fn1` and the estimation results (`fn2` and `fn3`).
The step widths for mcmc are determined as standard deviations / `error_scale`.

* `fn1`: Seafloor transponder position file (`fn1="pxp-ini.inp"` by default)
* `fn2`: NTD estimation result file obtained by `pos_array_all` (`fn2="residual.out"` by default)
* `fn3`: Estimation result file for all parameters obtained by `pos_array_all` (`fn3="solve.out"` by default)
* `error_scale`: Scaling factor for step widths (`error_scale=5.0` by default)
* `tscale`: Temporal scaling for time in the polynomial functions (time [sec] is converted into [hour]/`tscale`, `tscale=10` by default, this should correspond to the whole observational period in hours)

# Example
    make_initial(error_scale=6.0)
"""
function make_initial(;fn1="pxp-ini.inp"::String,fn2="residual.out"::String,fn3="solve.out"::String,fno="initial.inp"::String,error_scale = 5.0,tscale=10.0)
  println(stderr," === Make initial for static_array_mcmcgrad  ===")
  # --- Start log
  time1 = now()
  place = pwd()
  # --- Read data
  println(stderr," --- Read files")
  numk, px, py, pz = read_pxppos(fn1)
  ts, nk, to, tc, tr = read_ntd(fn2)
  dat = DelimitedFiles.readdlm(fn3)
  num0 = size(dat)[1]
  a0 = dat[:,1]; e0 = dat[:,2]/error_scale
  # --- Processing
  t0 = findmin(ts)[1] # start time in observational hour
  dep = mean(abs.(pz)) / 1000.0 # mean depth
  println(stderr," t0: $t0, dep: $dep")
  tt = (ts .- t0)/60/60/tscale  # observational hours
  H = hcat(tt*0 .+1.0,tt,tt.^2,tt.^3,tt.^4)
  d = copy(to)
  # --- Inversion
  println(stderr," --- Inversion for long-term NTD")
  dcal, dres, a, e = simple_inversion(d,H)
  num = size(a)[1]
  e = e/error_scale
  # --- Make initial.inp
  open(fno,"w") do out
    # Initial.inp: initial_value mimum_value maximum_value step_width
    println(out,"$(a0[1]) -20.0 20.0 $(e0[1]) EW_disp.") # EW disp.
    println(out,"$(a0[2]) -20.0 20.0 $(e0[2]) NS_disp.") # NS disp.
    println(out,"$(a0[3]) -20.0 20.0 $(e0[3]) UD_disp.") # UD disp.
    println(out,"0.0 -20.0 20.0 1.e-6 S-Gradient_EW") # EW shallow gradient
    println(out,"0.0 -20.0 20.0 1.e-6 S-Gradient_NS") # NS shallow gradient
    println(out,"0.65 0 $dep 0.005 Gradient_depth") # Gradient depth
    for i in 1:num
      println(out,"$(a[i]) -20.0 20.0 $(e[i]) L-NTD_$i") # Long-term NTD by 4d polynomial functions
    end
    println(out,"-4.0 -20.0 20.0 0.005 Scale_1") # Scaling factor for MCMC-1
    println(out,"-3.5 -20.0 20.0 0.005 Scale_2") # Scaling factor for MCMC-2
    for i in 4:num0
      ii = i - 3
      println(out,"$(a0[i]) -20.0 20.0 $(e0[i]) S-NTD_$ii") # 3d B-spline functions for detailed NTD
    end
  end
    
# --- Close process --- #
  time2 = now()
  println(stderr," Start time:",time1)
  println(stderr," Finish time:",time2)
  println(stderr," Check $fno")
end
