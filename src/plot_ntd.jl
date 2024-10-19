#using Plots
#using DelimitedFiles

export plot_ntd
"""
    plot_ntd(ntdrange,resrange; fno,fn,autoscale,plot_size,lmargin,rmargin,tmargin,bmargin,show,ms,bmargin0)

Make a figure of the time-series of the projected travel-time residuals in the nadir direction and the travel-time residuals removing the modeled NTD.

* `ntdrange`: Range of Y-axis for the projected travel-time residuals in the nadir direction
* `resrange`: Range of Y-axis for the travel-time residuals removing the modeled NTD
* `fn`: Input file name (`fn="residual.out"` by default)
* `fno`: Output figure file name (`fno="ntd.pdf"`)
* `plot_size`: Figure size (`plot_size=(650,650)` by default)
* `autoscale`: If `autoscale=true` (default), the plot range is automatically determined. If `autoscale=false`, the plot range of Y-axis is fixed by `ntdrange` and `resrange`
* `lmargin`: Plot margin for the left edge (`lmargin=2.5` by default)
* `rmargin`: Plot margin for the right edge (`rmargin=1.0` by default)
* `tmargin`: Plot margin for the top edge (`tmargin=1.0` by default)
* `bmargin`: Plot margin for the bottom edge (`bmargin=1.0` by default)
* `bmargin0`: Plot margin for the bottom edges of upper two panels (`bmargin0=-4.0` by default)
* `ms`: Plotted marker size (`ms=4` by default)

# Example
    plot_ntd(fno="ntd.pdf")
"""
function plot_ntd(ntdrange=(-3,3),resrange=(-1,1); fno="ntd.pdf"::String,fn="residual.out"::String,autoscale=true,plot_size=(650,650),lmargin=2.5, rmargin=1.0, tmargin=1.0, bmargin=1.0, show=false,ms=4::Int64,bmargin0=-4)
  dat0 = DelimitedFiles.readdlm(fn)
  num = size(dat0)[1]
  dat = unixsort2(dat0,1,2)
  t0 = dat[1,1]
  numk = round(Int64,findmax(dat[1:num,2])[1])
  if autoscale == true
    p0 = plot(ylabel="P-TT residual [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin0),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin),label="",xformatter=_->"")
    p1 = plot(xlabel="Time [hour]",ylabel="P-TT residual - NTD [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, rmargin),legend=:none)
  else
    p0 = plot(ylim=ntdrange,ylabel="P-TT residual [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin0),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin),label="",xformatter=_->"")
    p1 = plot(ylim=resrange,xlabel="Time [hour]",ylabel="P-TT residual -NTD [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, rmargin),legend=:none)
  end
  for k in 1:numk
    t = Vector{Float64}(undef, 0)
    d = Vector{Float64}(undef, 0)
    dc = Vector{Float64}(undef, 0)
    dr = Vector{Float64}(undef, 0)
    i0 = 1
    for i in 1:num
      if round(Int64,dat[i,2]) == k
        push!(t,((dat[i,1] - t0) / (60*60)))
        push!(d,(dat[i,3] * 1000))
        push!(dc,(dat[i,4] * 1000))
        push!(dr,(dat[i,5] * 1000))
        i0 += 1
      end
    end
    scatter!(p0,t,d,markershape=:cross,label="Transponder $k",markersize=ms)
    scatter!(p1,t,dr,markershape=:cross,markersize=ms)
  end
  t = (dat[1:num,1] .- t0) / (60*60)
  dc = dat[1:num,4] * 1000
  scatter!(p0,t,dc,markerstrokecolor=:black,mc=:black,markershape=:circle,markerstrokewidth=0,markersize=1,label="NTD")
  plts = plot(p0,p1,layout=(2,1), size=plot_size,framestyle=:box)
  if show == false
    savefig(plts,fno)
  else
    gui(plts)
  end
end
