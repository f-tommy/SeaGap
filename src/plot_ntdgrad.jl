#using Plots
#using DelimitedFiles

export plot_ntdgrad
"""
    plot_ntdgrad(ntdrange,resrange; fno,fn,autoscale,plot_size,lmargin,rmargin,tmargin,bmargin,show,ms,bmargin0)

Make a figure of the time-series of the projected travel-time residuals in the nadir direction, those without the deep gradients, and the travel-time residuals removing the modeled NTD.

* `ntdrange`: Range of Y-axis for the projected travel-time residuals in the nadir direction
* `resrange`: Range of Y-axis for the travel-time residuals removing the modeled NTD
* `fn`: Input file name (`fn="residual_grad.out"` by default)
* `fno`: Output figure file name (`fno="ntdgrad.pdf"`)
* `plot_size`: Figure size (`plot_size=(650,850)` by default)
* `autoscale`: If `autoscale=true` (default), the plot range is automatically determined. If `autoscale=false`, the plot range of Y-axis is fixed by `ntdrange` and `resrange`
* `lmargin`: Plot margin for the left edge (`lmargin=2.5` by default)
* `rmargin`: Plot margin for the right edge (`rmargin=1.0` by default)
* `tmargin`: Plot margin for the top edge (`tmargin=1.0` by default)
* `bmargin`: Plot margin for the bottom edge (`bmargin=1.0` by default)
* `bmargin0`: Plot margin for the bottom edges of upper two panels (`bmargin0=-4.0` by default)
* `ms`: Plotted marker size (`ms=4` by default)

# Example
    plot_ntdgrad()
"""
function plot_ntdgrad(ntdrange=(-3,3),resrange=(-1,1); fno="ntdgrad.pdf"::String,fn="residual_grad.out"::String,autoscale=true,plot_size=(650,850),lmargin=2.5, rmargin=1.0, tmargin=1.0, bmargin=1.0, show=false,ms=4::Int64,bmargin0=-4)
  dat0 = DelimitedFiles.readdlm(fn)
  num = size(dat0)[1]
  dat = unixsort2(dat0,1,2)
  t0 = dat[1,1]
  numk = round(Int64,findmax(dat[1:num,2])[1])
  if autoscale == true
    p0 = plot(ylabel="NTD [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin0),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin),label="",xformatter=_->"")
    p1 = plot(ylabel="NTD [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin0),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin),label="",xformatter=_->"")
    p2 = plot(xlabel="Time [hour]",ylabel="Residual [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, rmargin),legend=:none)
  else
    p0 = plot(ylim=ntdrange,ylabel="NTD [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin0),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin),label="",xformatter=_->"")
    p1 = plot(ylim=ntdrange,ylabel="NTD [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin0),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin),label="",xformatter=_->"")
    p2 = plot(ylim=resrange,xlabel="Time [hour]",ylabel="Residual [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, rmargin),legend=:none)
  end
  for k in 1:numk
    t = Vector{Float64}(undef, 0)
    d0 = Vector{Float64}(undef, 0)
    d = Vector{Float64}(undef, 0)
    dr = Vector{Float64}(undef, 0)
    i0 = 1
    for i in 1:num
      if round(Int64,dat[i,2]) == k
        push!(t,((dat[i,1] - t0) / (60*60)))
        push!(d0,(dat[i,6] * 1000))
        push!(d,((dat[i,6]-dat[i,8]) * 1000))
        push!(dr,(dat[i,13] * 1000))
        i0 += 1
      end
    end
    scatter!(p0,t,d0,ylabel="NTD [msec]",markershape=:cross,label="Transponder $k",markersize=ms)
    scatter!(p1,t,d,ylabel="NTD [msec]",markershape=:cross,label=:none,markersize=ms)
    scatter!(p2,t,dr,xlabel="Time [hour]",ylabel="Residual [msec]",markershape=:cross,label="$k",markersize=ms)
  end
  t = (dat[:,1] .- t0) / (60*60)
  d0 = dat[:,7] * 1000
  d1 = dat[:,9] * 1000
  d2 = dat[:,12] * 1000
  dl = dat[:,10] * 1000
  scatter!(p0,t,d1,mc=:black,markershape=:circle,markerstrokewidth=0,markersize=1,label="S-NTD & deep gradient")
  plot!(p1,t,d0,lc=:black,label="S-NTD")
  plot!(p1,t,dl,lc=:red,label="L-NTD")
  plot!(p1,t,d2,lc=:blue,label="L-NTD & shallow gradient")
  plts = plot(p0,p1,p2,layout=(3,1), size=plot_size,framestyle=:box)
  if show == false
    savefig(plts,fno)
  else
    gui(plts)
  end
end
