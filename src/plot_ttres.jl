#using Plots
#using DelimitedFiles

export plot_ttres
"""
    plot_ttres(resrange; autoscale,fn,fn0,plot_size,lmargin,rmargin,tmargin,bmargin,show,ms)

Make a figure of travel-time residuals obtained by `ttres()`

* `resrange`: Plot range of Y-sxis
* `autoscale`: If `autoscale=true` (default), the plot range is automatically determined. If `autoscale=false`, the plot range of Y-axis is fixed by `resrange`.
* `fn`: Input file name (`fn="ttres.out"` in default)
* `fno`: Output figure name (`fno="ttres.pdf"` in default)
* `plot_size`: Figure size (`plot_size=(650,1200)` in default)
* `lmargin`: Plot margin for the left edge (`lmargin=6.0` in default)
* `rmargin`: Plot margin for the right edge (`rmargin=1.0` in default)
* `tmargin`: Plot margin for the top edge (`tmargin=1.0` in default)
* `bmargin`: Plot margin for the bottom edge (`bmargin=1.0` in default)
* `show`: if `show=true`, a figure is temporally shown; if false, the figure is saved as `fno` (`show=false` in default)
* `ms`: Plotted marker size (`ms=7` in default)

# Example

"""
function plot_ttres(resrange=(-3,3); autoscale=true::Bool,fn="ttres.pdf"::String,fn0="ttres.out"::String,plot_size=(650,1200),lmargin=6.0, rmargin=1.0, tmargin=1.0, bmargin=1.0, show=false::Bool,ms=7::Int64)
  dat0 = DelimitedFiles.readdlm(fn0)
  num = size(dat0)[1]
  dat = unixsort2(dat0,3,2)
  t0 = dat[1,3]
  numk = round(Int64,findmax(dat[1:num,2])[1])
  p = [] 
  if numk == 1
    if autoscale == true
      push!(p,plot(xlabel="Time [hours]", ylabel="Residual [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin)))
    else
      push!(p,plot(ylim=resrange, xlabel="Time [hours]", ylabel="Residual [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin)))
    end
  else
    for k in 1:numk
      if k == 1
        if autoscale == true
          push!(p,plot(ylabel="Residual [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, 0),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin)))
        else
          push!(p,plot(ylim=resrange,ylabel="Residual [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, 0),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin)))
        end
      elseif k == numk
        if autoscale == true
          push!(p,plot(xlabel="Time [hours]",ylabel="Residual [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, rmargin)))
        else
          push!(p,plot(ylim=resrange,xlabel="Time [hours]",ylabel="Residual [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, rmargin)))
        end
      else
        if autoscale == true
          push!(p,plot(ylabel="Residual [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, 0),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, rmargin)))
        else
          push!(p,plot(ylim=resrange,ylabel="Residual [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, 0),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, rmargin)))
        end
      end
    end
  end
  for k in 1:numk
    t = Vector{Float64}(undef, 0)
    dr = Vector{Float64}(undef, 0)
    i0 = 1
    for i in 1:num
      if round(Int64,dat[i,2]) == k
        push!(t,((dat[i,3] - t0) / (60*60)))
        push!(dr,(dat[i,7] * 1000))
        i0 += 1
      end
    end
    scatter!(p[k],t,dr,markershape=:cross,label="$k",markersize=ms)
  end
  plts = plot(p...,layout=(numk,1), size=plot_size)
  if show == false
    savefig(plts,fn)
  else
    gui(plts)
  end
end
