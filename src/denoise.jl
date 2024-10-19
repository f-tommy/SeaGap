#using Plots
#using DelimitedFiles
#using Statistics
#using Printf

export plot_denoise, denoise
"""
    plot_denoise(;resrange,resrange2,autoscale,fn,fno,plot_size,lmargin,rmargin,tmargin,bmargin,show,ms)

Make a figure of `fno` from the `denoise()` output file.

* `resrange`: Vertical range for the vertically-projected travel-time residual [ms] (with NTD) (`resrange=(-3,3)` by default)
* `resrange2`: Vertical range for the vertically-projected travel-time residual [ms] (excluding NTD) (`resrange2=(-1,1)` by default)
* `autoscale`: if `autoscale=true`, the vertical ranges are automatically determined; if `autoscale=false`, the vertical ranges are fixed to `resrange` and `resrange2` (`autoscale=true` by default)
* `plot_size`: Figure size (`plot_size=(1200,1200)` by default)
* `fn`: Input file name (`fn="denoise.out"` by defualt)
* `fno`: Output file name (`fno="denoise.png"` by defualt)
* `lmargin`: Plot margin for the left edge (`lmargin=6.0` by default)
* `rmargin`: Plot margin for the right edge (`rmargin=1.0` by default)
* `tmargin`: Plot margin for the top edge (`tmargin=1.0` by default)
* `bmargin`: Plot margin for the bottom edge (`bmargin=1.0` by default)
* `show`: if `show=true`, a figure is temporally shown; if false, the figure is save as `fno` (`show=false` by default)
* `ms`: Plotted marker size (`ms=6` by default)
"""
function plot_denoise(;resrange=(-3,3),resrange2=(-1,1), autoscale=true::Bool,fno="denoise.png"::String,fn="denoise.out"::String, plot_size=(1200,1200),lmargin=6.0, rmargin=1.0, tmargin=1.0, bmargin=1.0, show=false::Bool, ms=6::Int64)
  dat0 = DelimitedFiles.readdlm(fn)
  num = size(dat0)[1]
  dat = unixsort2(dat0,3,2)
  t0 = dat[1,3]
  numk = round(Int64,findmax(dat[1:num,2])[1])
  p1 = [] 
  p2 = [] 
  println(stderr,"Set base-figure")
  if numk == 1
    if autoscale == true
      push!(p1,plot(xlabel="Time [hours]", ylabel="Residual [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, 0.0)))
      push!(p2,plot(xlabel="Time [hours]", left_margin=Plots.Measures.Length(:mm, 0.0),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin)))
    else
      push!(p1,plot(ylim=resrange, xlabel="Time [hours]", ylabel="Residual [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, 0.0)))
      push!(p2,plot(ylim=resrange2, xlabel="Time [hours]", left_margin=Plots.Measures.Length(:mm, 0.0),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin)))
    end
  else
    for k in 1:numk
      if k == 1
        if autoscale == true
          push!(p1,plot(ylabel="Residual [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, 0),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, 0.0)))
          push!(p2,plot(left_margin=Plots.Measures.Length(:mm, 0.0),bottom_margin=Plots.Measures.Length(:mm, 0),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin)))
        else
          push!(p1,plot(ylim=resrange,ylabel="Residual [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, 0),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, 0.0)))
          push!(p2,plot(ylim=resrange2,left_margin=Plots.Measures.Length(:mm, 0.0),bottom_margin=Plots.Measures.Length(:mm, 0),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin)))
        end
      elseif k == numk
        if autoscale == true
          push!(p1,plot(xlabel="Time [hours]",ylabel="Residual [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, 0.0)))
          push!(p2,plot(xlabel="Time [hours]",left_margin=Plots.Measures.Length(:mm, 0.0),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, rmargin)))
        else
          push!(p1,plot(ylim=resrange,xlabel="Time [hours]",ylabel="Residual [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, 0.0)))
          push!(p2,plot(ylim=resrange2,xlabel="Time [hours]",left_margin=Plots.Measures.Length(:mm, 0.0),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, rmargin)))
        end
      else
        if autoscale == true
          push!(p1,plot(ylabel="Residual [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, 0),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, 0)))
          push!(p2,plot(left_margin=Plots.Measures.Length(:mm, 0),bottom_margin=Plots.Measures.Length(:mm, 0),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, rmargin)))
        else
          push!(p1,plot(ylim=resrange,ylabel="Residual [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, 0),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, 0)))
          push!(p2,plot(ylim=resrange2,left_margin=Plots.Measures.Length(:mm, 0),bottom_margin=Plots.Measures.Length(:mm, 0),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, rmargin)))
        end
      end
    end
  end
  println("Plot residuals")
  for k in 1:numk
    for j in 0:1
      t = Vector{Float64}(undef, 0)
      r = Vector{Float64}(undef, 0)
      dr = Vector{Float64}(undef, 0)
      for i in 1:num
        if round(Int64,dat[i,2]) == k
          if round(Int64,dat[i,7]) == j
            push!(t,((dat[i,3] - t0) / (60*60)))
            push!(r,(dat[i,4] * 1000))
            push!(dr,(dat[i,6] * 1000))
          end
        end
      end
      if j == 1
        scatter!(p1[k],t,r,markershape=:cross,label=:none, mc=:red, markersize=ms)
        scatter!(p2[k],t,dr,markershape=:cross,label=:none, mc=:red, markersize=ms)
      else
        scatter!(p1[k],t,r,markershape=:cross,label="$k", mc=:blue, markersize=ms)
        scatter!(p2[k],t,dr,markershape=:cross,label="$k", mc=:blue, markersize=ms)
      end
    end
    t = Vector{Float64}(undef, 0)
    sr = Vector{Float64}(undef, 0)
    for i in 1:num
      if round(Int64,dat[i,2]) == k
        push!(t,((dat[i,3] - t0) / (60*60)))
        push!(sr,(dat[i,5] * 1000))
      end
    end
    plot!(p1[k],t,sr,label=:none, lc=:magenta, lw=2)
  end
  p3 = []
  for k in 1:numk
    push!(p3,plot(p1[k],p2[k],layout=(1,2)))
  end
  plts = plot(p3...,layout=(numk,1), size=plot_size)
  if show == false
    savefig(plts,fno)
  else
    gui(plts)
  end
end

"""
    denoise(lat,TR_DEPTH,resrange,resrange2; method,autoscale,k,n,sigma,save,prompt,fn1,fn2,fn3,fn4,fn0,plot_size,lmargin,rmargin,tmargin,bmargin,show,fno1,fno2)

Calculate travel-time residual, estimate smoothed travel-time residuals by `n` running `method` filter, exclude outliers beyond `sigma` Std, and plot them by `plot_denoise()`. The denoised observational file is rewritten in `fn4`.

* `lat`: Site latitude
* `TR_DEPTH`: Transducer depth from the sea-surface
* `method`: Method of running filter ("mean" or "method"; `method="median"` by default)
* `k`: if `k=0`, `denoise()` is performed for all transponders; if `k` >= 1, `denoise()` is performed for the `k`th transponder. 
* `n`: Window size for the running filter
* `sigma`: Outlier threshold
* `save`: if `save=true`, the original observation data file `fn4` is renamed and saved as `fn0` (`save=true` by default)
* `prompt`: if `prompt=true`, confirmation message is shown; if false, the denoised observation file is forcely saved (`prompt=true` by default)  

* `fn1`: GNSS antenna-transducer offset (`fn1="tr-ant.inp"` by default)
* `fn2`: Initial transponders position (`fn2="pxp-ini.inp"` by default)
* `fn3`: Initial sound speed structure (`fn3="ss_prof.inp"` by default)
* `fn4`: Observational file (`fn4="obsdata.inp"` by default)
* `fn0`: if `save=true`, `fn4` is saved (`fn0="obsdata.inp_org"` by default)

* `resrange`: Vertical range for the vertically-projected travel-time residual [ms] (with NTD) (`resrange=(-3,3)` by default)
* `resrange2`: Vertical range for the vertically-projected travel-time residual [ms] (excluding NTD) (`resrange2=(-1,1)` by default)
* `autoscale`: if `autoscale=true`, the vertical ranges are automatically determined; if `autoscale=false`, the vertical ranges are fixed to `resrange` and `resrange2` (`autoscale=true` by default)
* `plot_size`: Figure size (`plot_size=(1200,1200)` by default)
* `lmargin`: Plot margin for the left edge (`lmargin=6.0` by default)
* `rmargin`: Plot margin for the right edge (`rmargin=1.0` by default)
* `tmargin`: Plot margin for the top edge (`tmargin=1.0` by default)
* `bmargin`: Plot margin for the bottom edge (`bmargin=1.0` by default)
* `show`: if `show=true`, a figure is temporally shown; if false, the figure is save as `fno` (`show=true` by default)
* `ms`: Plotted marker size (`ms=6` by default)

* `fno1`: Output text file for travel-time residuals
* `fno2`: Output figure name

# Example
    denoise(lat,TR_DEPTH,k=0,n=7,sigma=4.0,method="median")
"""
function denoise(lat,TR_DEPTH::Vector{Float64},resrange=(-3,3),resrange2=(-1,1); method="median"::String,autoscale=true::Bool,k=0,n=15,sigma=4.0,save=true::Bool,prompt=true::Bool,fn1="tr-ant.inp"::String, fn2="pxp-ini.inp"::String, fn3="ss_prof.inp"::String, fn4="obsdata.inp"::String, fn0="obsdata.inp_org"::String, plot_size=(1200,1200),lmargin=6.0, rmargin=1.0, tmargin=1.0, bmargin=1.0, show=true::Bool, fno1="denoise.out"::String, fno2="denoise.png"::String,ms=6::Int64)
  # Calcualte travel-time residuals
  nv,kv,t1,t2,tp,tc,tr,vert = ttres(lat,TR_DEPTH, fn1=fn1, fn2=fn2, fn3=fn3, fn4=fn4)
  dat = hcat(nv,kv,t1,t2,tp,tc,tr,vert)
  # Set array and parameters
  odat = DelimitedFiles.readdlm(fn4)
  println(stderr,"Set array and parameters for denoise")
  num = size(dat)[1]
  dat0 = unixsort2(dat,3,2)
  t0 = dat0[1,3]
  numk = round(Int64,findmax(dat[1:num,2])[1])
  if k > numk
    error("denoise.jl: k is incorrectly specified.")
  end
  if save ==true 
    cp("$fn4","$fn0",force=true)
    println(stderr,"Copy $fn4 -> $fn0")
  end
  # Running filter
  id = zeros(Int64,num)
  open(fno1,"w") do out
  for j in 1:numk
    println(stderr,"Running filter for Transponder $j")
    i0 = Vector{Int64}(undef, 0)
    t = Vector{Float64}(undef, 0)
    r = Vector{Float64}(undef, 0)
    for i in 1:num
      if round(Int64,dat[i,2]) == j
        push!(i0,round(Int64,dat[i,1]))
        push!(t,(dat[i,3]))
        push!(r,(dat[i,7]))
      end
    end
    numi = length(r)
    println(stderr,"  $numi")
    if method == "median"
      sr = runmed(r,n)
    else
      sr = runave(r,n)
    end
    dr = r - sr
    sd = std(dr)
    for i in 1:numi
      id0 = 0
      if j == k || k == 0
        if abs(dr[i]) > sigma*sd
          ii = i0[i]
          id[ii] = 1
          id0 = 1
        end
      end
      println(out,"$(i0[i]) $j $(t[i]) $(r[i]) $(sr[i]) $(dr[i]) $id0")
    end
  end
  end
  # Make a new observational data file
  println(stderr,"Make denoised observational data file: $fn4")
  open("tmp","w") do out
  for i in 1:num
    if id[i] == 0
      @printf(out,"%i %2.6f %10.6f %6.6f %6.6f %6.6f %5.5f %3.5f %3.5f %10.6f %6.6f %6.6f %6.6f %5.5f %3.5f %3.5f %d %d\n",odat[i,1],odat[i,2],odat[i,3],odat[i,4],odat[i,5],odat[i,6],odat[i,7],odat[i,8],odat[i,9],odat[i,10],odat[i,11],odat[i,12],odat[i,13],odat[i,14],odat[i,15],odat[i,16],odat[i,17],odat[i,18])
    end
  end
  end
  # Plot
  println(stderr,"Plot denoised time-series")
  plot_denoise(fno=fno2,fn=fno1, resrange=resrange,resrange2=resrange2,plot_size=plot_size,lmargin=lmargin, rmargin=rmargin, tmargin=tmargin, bmargin=bmargin, show=show,ms=ms) 
  if prompt == true
    q = Base.prompt("Do you accept the denoise processing? (yes/no)")
    if q == "yes"
      mv("tmp",fn4,force=true)
    else
      rm("tmp")
    end
  else
     mv("tmp",fn4,force=true)
  end
end 
