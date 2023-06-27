#using Plots
#using DelimitedFiles

export plot_track, plot_timetrack
"""
    plot_track(xrange,yrange;autoscale,fn1,fn2,fno,plot_size,lmargin,tmargin,bmargin,rmargin,show,ms1,ms2,gfs,cb)

* `xrange`: Range of EW component in meters when `autoscale=false`
* `yrange`: Range of NS component in meters when `autoscale=false`
* `autoscale`: If `autoscale=true`, the plot range is automatically defined
* `fn1`: File name of the seafloor transponder positions (`fn1="pxp-ini.inp"` by default)
* `fn2`: Travel-time residual file name (`fn2="obsdata.inp"` by default)
* `fno`: Output figure name (`fno="track.pdf"` by default)
* `show`: if `show=true`, a figure is shown on REPL and is not saved as a file (`show=false` by default)
* `plot_size`: Figure size (`plot_size=(600,500)` by default)
* `lmargin`: Plot margin for the left edge (`lmargin=2.5` by default)
* `rmargin`: Plot margin for the right edge (`rmargin=1.0` by default)
* `tmargin`: Plot margin for the top edge (`tmargin=1.0` by default)
* `bmargin`: Plot margin for the bottom edge (`bmargin=1.0` by default)
* `ms1`: Plotted marker size for transponders (`ms1=10` by default)
* `ms2`: Plotted marker size for track (`ms2=5` by default)
* `gfs`: Label fontsize (`gfs=12` by default)
* `cb`: Color scale (`cb=:lightrainbow` by default)

# Example
    plot_track((-3000,3000),(-3000,3000),autoscale=false,fn1="pxp-ini.inp",fn2="obsdata.inp",fno="track.pdf")
"""
function plot_track(xrange=(-3000,3000),yrange=(-3000,3000);autoscale=true::Bool,fn1="pxp-ini.inp"::String,fn2="obsdata.inp"::String,fno="track.pdf"::String,plot_size=(600,500),lmargin=2.5,tmargin=1.0,bmargin=1.0,rmargin=1.0,show=false::Bool,ms1=10::Int64,ms2=5::Int64,gfs=12::Int64,cb=:rainbow)
  a = DelimitedFiles.readdlm(fn1)
  numk = size(a)[1]
  px = a[1:numk,1]
  py = a[1:numk,2]
  a = DelimitedFiles.readdlm(fn2)
  num = size(a)[1]
  x = a[1:num,4]
  y = a[1:num,5]
  ts = a[1,3]
  te = a[num,3]
  t = (a[1:num,3].-ts)/60/60
  dt = floor(Int64,(te-ts)/4/60/60)
  te = (te-ts)/60/60
  ts = 0.0
  cpt = cgrad(cb,[ts,dt,te])
  p = scatter(px,py,markershape=:utriangle,markersize=ms1,c=RGB(0.8,0.8,0.8))
  if autoscale == true
    scatter!(p,x,y,marker_z=t,markershape=:cross,c=cpt,xlabel="Easting [m]",ylabel="Northing [m]",aspect_ratio = 1,legend=false,colorbar=true, framestyle=:box,colorbar_title="Time [hour]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin),markersize=ms2,guidefontsize=gfs)
  else
    scatter!(p,x,y,marker_z=t,markershape=:cross,c=cpt,xlabel="Easting [m]",ylabel="Northing [m]",xlims=xrange,ylims=yrange,aspect_ratio = 1,legend=false,colorbar=true,colorbar_title="Time [hour]",framestyle=:box,size=plot_size,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin),markersize=ms2,guidefontsize=gfs)
  end
  if show == false
    savefig(p,fno)
  else
    gui(p)
  end
end

"""
    plot_timetrack(xrange=(-3000,3000),yrange=(-3000,3000),zrange=(0,10) ;autoscale=true::Bool,fn="obsdata.inp"::String,fno="timetrack.pdf"::String,plot_size=(600,800),lmargin=4.0,tmargin=1.0,bmargin=0.5,rmargin=1.0,bmargin0=-3.0,show=false::Bool,ms=5::Int64,gfs=12::Int64)

Make a figure on time-series of the sea-surface platform positions

* `xrange`, `yrange`, and `zrange`: 
* `fn`: Input observational file (`fn="obsdata.inp"` by default)
* `fno`: Output figure name (`fno="timetrack.pdf"`)
* `autoscale`: If `autoscale=true` (default), the plot range is automatically determined. If `autoscale=false`, the plot range of Y-axis is fixed by `xrange`, `yrange`, and `zrange`.
* `plot_size`: Figure size (`plot_size=(600,800)` by default)
* `lmargin`: Plot margin for the left edge (`lmargin=4.0` by default)
* `rmargin`: Plot margin for the right edge (`rmargin=1.0` by default)
* `tmargin`: Plot margin for the top edge (`tmargin=1.0` by default)
* `bmargin`: Plot margin for the bottom edge (`bmargin=0.5` by default)
* `bmargin0`: Plot margin for the bottom edges of upper two panels (`bmargin0=-3.0` by default)
* `show`: if `show=true`, a figure is temporally shown; if false, the figure is saved as `fno` (`show=false` by default)
* `ms`: Plotted marker size (`ms=5` by default)
* `gfs`: Fontsize for label (gudefontsize: `gfs=12` by default)

# Example
    plot_timetrack(fno="timetrack.png")
"""
function plot_timetrack(xrange=(-3000,3000),yrange=(-3000,3000),zrange=(0,10) ;autoscale=true::Bool,fn="obsdata.inp"::String,fno="timetrack.pdf"::String,plot_size=(600,800),lmargin=4.0,tmargin=1.0,bmargin=0.5,rmargin=1.0,bmargin0=-3.0,show=false::Bool,ms=5::Int64,gfs=12::Int64)
  a = DelimitedFiles.readdlm(fn)
  num = size(a)[1]
  x = a[1:num,4]
  y = a[1:num,5]
  z = a[1:num,6]
  ts = a[1,3]
  te = a[num,3]
  t = (a[1:num,3].-ts)/60/60
  te = (te-ts)/60/60
  ts = 0.0
  if autoscale == false
    p1 = scatter(t,x,xlims=(ts,te),ylims=xrange,ylabel="Easting [m]",framestyle=:box,legend=:none,markershape=:cross,markersize=ms,label="",xformatter=_->"",top_margin=Plots.Measures.Length(:mm,tmargin),bottom_margin=Plots.Measures.Length(:mm,bmargin0))
    p2 = scatter(t,y,xlims=(ts,te),ylims=yrange,ylabel="Northing [m]",framestyle=:box,legend=:none,markershape=:cross,markersize=ms,label="",xformatter=_->"",top_margin=Plots.Measures.Length(:mm,0),bottom_margin=Plots.Measures.Length(:mm,bmargin0))
    p3 = scatter(t,z,xlims=(ts,te),ylims=zrange,ylabel="Height [m]",framestyle=:box,legend=:none,markershape=:cross,markersize=ms,top_margin=Plots.Measures.Length(:mm,0),bottom_margin=Plots.Measures.Length(:mm,bmargin))
  else
    p1 = scatter(t,x,xlims=(ts,te),ylabel="Easting [m]",framestyle=:box,legend=:none,markershape=:cross,markersize=ms,label="",xformatter=_->"",top_margin=Plots.Measures.Length(:mm,0),bottom_margin=Plots.Measures.Length(:mm,bmargin0))
    p2 = scatter(t,y,xlims=(ts,te),ylabel="Northing [m]",framestyle=:box,legend=:none,markershape=:cross,markersize=ms,label="",xformatter=_->"",top_margin=Plots.Measures.Length(:mm,0),bottom_margin=Plots.Measures.Length(:mm,bmargin0))
    p3 = scatter(t,z,xlims=(ts,te),xlabel="Time [hour]",framestyle=:box,ylabel="Height [m]",legend=:none,markershape=:cross,markersize=ms,top_margin=Plots.Measures.Length(:mm,0),bottom_margin=Plots.Measures.Length(:mm,bmargin))
  end
  plt = plot(p1,p2,p3,layout=(3,1),size=plot_size,guidefontsize=gfs,left_margin=Plots.Measures.Length(:mm,lmargin),right_margin=Plots.Measures.Length(:mm, rmargin))
  if show == false
    savefig(plt,fno)
  else
    gui(plt)
  end
end
