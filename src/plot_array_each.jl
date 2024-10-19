#using Plots
#using DelimitedFiles

export plot_time_array_each, plot_map_array_each

"""
    plot_map_array_each(xrange,yrange; autoscale,fn,fno,plot_size,lmargin,tmargin,bmargin,rmargin,show,ms,gfs,col_num)

Make a figure plotting the estimated array displacements obtained by `pos_array_each()` in a horizontal map.

* `xrange` and `yrange`: EW and NS ranges for plot [m]
* `autoscale`: If `autoscale=true` (default), the plot range is automatically determined. If `autoscale=false`, the plot range is fixed by `xrange` and `yrange`.
* `fn`: Input file name (`fn="array_each.out"` by default)
* `fno`: Output figure name (`fno="map_array_each.pdf"` by default)
* `plot_size`: Figure size (`plot_size=(600,500)` by default)
* `lmargin`: Plot margin for the left edge (`lmargin=2.5` by default)
* `rmargin`: Plot margin for the right edge (`rmargin=1.0` by default)
* `tmargin`: Plot margin for the top edge (`tmargin=1.0` by default)
* `bmargin`: Plot margin for the bottom edge (`bmargin=1.0` by default)
* `show`: if `show=true`, a figure is temporally shown; if false, the figure is save as `fno` (`show=false` by default)
* `ms`: Plotted marker size (`ms=5` by default)
* `gfs`: Fontsize for label (gudefontsize: `gfs=12` by default)
* `col_num`: If `col_num=1` (default), the plot is colored by the observation time. If `col_num=2`, the plot is colored by number of the used observational data. If `col_num=0`, the plot is colored by blue.

# Example
    SeaGap.plot_map_array_each((-0.5,0.5),(-0.5,0.5),autoscale=false)

"""
function plot_map_array_each(xrange=(-1,1),yrange=(-1,1);autoscale=true::Bool,fn="array_each.out"::String,fno="map_array_each.pdf"::String,plot_size=(600,500),lmargin=2.5,tmargin=1.0,bmargin=1.0,rmargin=1.0,show=false::Bool,ms=5::Int64,gfs=12::Int64,col_num=1::Int64)
  a = DelimitedFiles.readdlm(fn)
  x = a[:,3]
  y = a[:,4]
  if col_num == 1
    z = (a[:,1].-a[1,1])/60/60
    ct = "Time [hour]"
  elseif col_num > 1
    z = a[:,col_num]
    if col_num == 2
      ct = "Number of data"
    else
      ct = ""
    end
  end
  if autoscale == true
    if col_num == 0
      p = scatter(x,y,markershape=:cross,xlabel="Easting [m]",ylabel="Northing [m]",aspect_ratio = 1,legend=false,c=:blue,framestyle=:box,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin),markersize=ms,guidefontsize=gfs)
    else
      p = scatter(x,y,zcolor=z,markershape=:cross,c=:rainbow,xlabel="Easting [m]",ylabel="Northing [m]",aspect_ratio = 1,legend=false,colorbar=true, framestyle=:box,colorbar_title=ct,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin),markersize=ms,guidefontsize=gfs)
    end
  else
    if col_num == 0
      p = scatter(x,y,markershape=:cross,xlabel="Easting [m]",ylabel="Northing [m]",xlims=xrange,ylims=yrange,aspect_ratio = 1,c=:blue,legend=false,framestyle=:box,size=plot_size,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin),markersize=ms,guidefontsize=gfs)
    else
      p = scatter(x,y,zcolor=z,markershape=:cross,c=:rainbow,xlabel="Easting [m]",ylabel="Northing [m]",xlims=xrange,ylims=yrange,aspect_ratio = 1,legend=false,colorbar=true,colorbar_title=ct,framestyle=:box,size=plot_size,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin),markersize=ms,guidefontsize=gfs)
    end
  end
  if show == false
    savefig(p,fno)
  else
    gui(p)
  end
end

"""
    plot_time_array_each(EW_range,NS_range, ntdrange; autoscale,fn,fno,plot_size,lmargin,tmargin,bmargin,rmargin,bmargin0,show,ms,gfs)

Make a figure plotting the estimated array displacements obtained by `pos_array_each()` in time-series.

* `EW_range`, `NS_range`, and `ntdrange`: Plot range of Y-sxis for EW [m], NS [m], and NTD [ms] components
* `autoscale`: If `autoscale=true` (default), the plot range is automatically determined. If `autoscale=false`, the plot range of Y-axis is fixed by `EW_range`, `NS_range`, and `ntdrange`.
* `fn`: Input file name (`fn="array_each.out"` by default)
* `fno`: Output figure name (`fno="time_array_each.pdf"` by default)
* `plot_size`: Figure size (`plot_size=(600,800)` by default)
* `lmargin`: Plot margin for the left edge (`lmargin=3.5` by default)
* `rmargin`: Plot margin for the right edge (`rmargin=1.0` by default)
* `tmargin`: Plot margin for the top edge (`tmargin=1.0` by default)
* `bmargin`: Plot margin for the bottom edge (`bmargin=1.0` by default)
* `bmargin0`: Plot margin for the bottom edges of upper two panels (`bmargin0=-4.0` by default)
* `show`: if `show=true`, a figure is temporally shown; if false, the figure is saved as `fno` (`show=false` by default)
* `ms`: Plotted marker size (`ms=6` by default)
* `gfs`: Fontsize for label (gudefontsize: `gfs=12` by default)

# Example
    SeaGap.plot_time_array_each(show=true)

"""
function plot_time_array_each(EW_range=(-1.5,1.5),NS_range=(-1.5,1.5),ntdrange=(-3,3);autoscale=true::Bool,fno="time_array_each.pdf"::String,fn="array_each.out"::String,plot_size=(600,800),lmargin=3.5,rmargin=1.0, tmargin=1.0, bmargin=1.0, bmargin0=-4.0,show=false,ms=6::Int64,gfs=12::Int64)
  a = DelimitedFiles.readdlm(fn)
  t0 = a[1,1]
  t = (a[:,1] .- t0) / (60*60)
  k = a[:,2]
  x = a[:,3]
  y = a[:,4]
  dt = a[:,5] * 1000
  if autoscale == true
    p0 = scatter(t,dt,ylabel="NTD [msec]",framestyle=:box,legend = :none,markershape=:cross,bottom_margin=Plots.Measures.Length(:mm,bmargin0),top_margin=Plots.Measures.Length(:mm, tmargin),markersize=ms,label="",xformatter=_->"")
    p1 = scatter(t,x,ylabel="Easting [m]",framestyle=:box,legend = :none,markershape=:cross,bottom_margin=Plots.Measures.Length(:mm,bmargin0),top_margin=Plots.Measures.Length(:mm, 0),markersize=ms,label="",xformatter=_->"")
    p2 = scatter(t,y,xlabel="Time [hour]",ylabel="Northing [m]",framestyle=:box,legend = :none,markershape=:cross,bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, 0),markersize=ms)
  else
    p0 = scatter(t,dt,ylabel="NTD [msec]",framestyle=:box,legend = :none,ylim=ntdrange,markershape=:cross,bottom_margin=Plots.Measures.Length(:mm,bmargin0),top_margin=Plots.Measures.Length(:mm, tmargin),markersize=ms,label="",xformatter=_->"")
    p1 = scatter(t,x,ylabel="Easting [m]",framestyle=:box,legend = :none,ylim=EW_range,markershape=:cross,bottom_margin=Plots.Measures.Length(:mm,bmargin0),top_margin=Plots.Measures.Length(:mm, 0),markersize=ms,label="",xformatter=_->"")
    p2 = scatter(t,y,xlabel="Time [hour]",ylabel="Northing [m]",framestyle=:box,legend = :none,ylim=NS_range,markershape=:cross,bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, 0),makersize=ms)
  end
  plts = plot(p0,p1,p2,layout=(3,1), size=plot_size, right_margin=Plots.Measures.Length(:mm, rmargin),left_margin=Plots.Measures.Length(:mm, lmargin),guidefontsize=gfs)
  if show == false
    savefig(plts,fno)
  else
    gui(plts)
  end
end


