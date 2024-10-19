#using Dates
#using Statistics
#using DelimitedFiles
#include("read_gnssa.jl")
#include("simple_inversion.jl")
# Usage :make_initial_gradv()

export make_initial_gradv
"""
    make_initial_gradv(NPB1,NPB2;fn1,fn2,fn3,fno,error_scale,gs,gd,sntd,lntd,scale,gdep0,sntdini,hp_err,hp_sntd,hp_gs,hp_gd,gsini)

Make an initial data file `fno` for `static_array_mcmc_all()` using initial transponder positions `fn1` and the estimation results (`fn2` and `fn3`).
The step widths for mcmc are determined as standard deviations / `error_scale`.

* `fn`: Estimation result file for all parameters obtained by `static_array_s` (`fn="solve.out"` by default)
* `error_scale`: Scaling factor for step widths (`error_scale=5.0` by default)
* `gs`: Step width for shallow gradients
* `gd`: Step width for gradient depth
* `sntd`: Step width for short-term NTDs. If `sntd` > 1, `stnd` is set following the standard deviations / `error_scale`.
* `lntd`: Step width for long-term NTDs. If `lntd` > 1, `ltnd` is set following the standard deviations / `error_scale`.
* `scale`: Step width for scaling factors for hyper-parameters
* `gdep0`: Initial value for gradient depth
* `hp_err`: Initial value for hyper-parameter of observational error
* `hp_sntd`: Initial value for hyper-parameter of S-NTD
* `hp_gs`: Initial value for hyper-parameter of shallow gradients
* `hp_gd`: Initial value for hyper-parameter of gradient depth
* `gsini`: if false, initial value of shallowe gradients are given as zero

# Example
    make_initial_gradv(5,100)
"""
function make_initial_gradv(NPB1=5::Int64,NPB2=100::Int64; fn="solve.out"::String,fno="initial.inp"::String,error_scale = 5.0, gs = 2.e-6, gd = 0.01, sntd = 1.e-6, lntd = 1.e-6, scale = 0.005, gdep0 = 0.0, sntdini = true::Bool, hp_err = -4.0, hp_sntd = -4.5, hp_gs = -4.5, hp_gd = -1.0, gsini = false::Bool,spc = false::Bool)
  println(stderr," === Make initial for static_array_mcmcgradv  ===")
  # --- Start log
  time1 = now()
  place = pwd()
  # --- Read data
  println(stderr," --- Read files")
  dat = DelimitedFiles.readdlm(fn)
  num0 = size(dat)[1]
  a0 = dat[:,1]; e0 = dat[:,2]/error_scale
  # --- Make initial.inp
  open(fno,"w") do out
    # Initial.inp: initial_value mimum_value maximum_value step_width
    println(out,"$(a0[1]) $(e0[1]) EW_disp.") # EW disp.
    println(out,"$(a0[2]) $(e0[2]) NS_disp.") # NS disp.
    println(out,"$(a0[3]) $(e0[3]) UD_disp.") # UD disp.
    if gsini == false
      println(out,"$(a0[4]) $gs S-Grad_EW") # EW shallow gradient
      println(out,"$(a0[5]) $gs S-Grad_NS") # NS shallow gradient
    else
      println(out,"0.0 $gs S-Grad_EW") # EW shallow gradient
      println(out,"0.0 $gs S-Grad_NS") # NS shallow gradient
    end
    println(out,"$gdep0 $gd G-Depth_EW") # EW gradient depth
    println(out,"$gdep0 $gd G-Depth_NS") # NS gradient depth
    println(out,"$hp_err $scale HP_error") # Scaling factor
    println(out,"$hp_sntd $scale HP_S-NTD") # B-spline constraint
    println(out,"$hp_gs $scale HP_S-Grad") # Norm constraint on shallow gradient
    println(out,"$hp_gd $scale HP_G-Depth") # Norm constraint on deep gradient
    for i in 1:NPB1
      ii = i + 5
      if lntd > 1.0
        println(out,"$(a0[ii]) $(e0[ii]) L-NTD_$i") # 3d B-spline functions for rough NTD
      else
        println(out,"$(a0[ii]) $lntd L-NTD_$i") # 3d B-spline functions for rough NTD
      end
    end
    if sntd > 1.0
      ee = mean(e0[6+NPB1:end])
      println(out,"0.0 $ee S-NTD") # 3d B-spline functions for detailed NTD
    else
      println(out,"0.0 $sntd S-NTD") # 3d B-spline functions for detailed NTD
    end
  end
    
# --- Close process --- #
  time2 = now()
  println(stderr," Start time:",time1)
  println(stderr," Finish time:",time2)
  println(stderr," Check $fno")
end
