#using Plots
#using DelimitedFiles

export plot_gradmap_gradv
"""
    plot_gradmap_gradv(xrange,yrange; autoscale,fn1,fn2,fno,plot_size,lmargin,tmargin,bmargin,rmargin,show,ms1,ms2,gfs,shape,cb)

Make a figure of sea-surface track with color of shallow gradient from the results of `static_array_mcmc_gradv`.

* `xrange`: Range of EW component in meters when `autoscale=false`
* `yrange`: Range of NS component in meters when `autoscale=false`
* `autoscale`: If `autoscale=true`, the plot range is automatically defined
* `fn1`: File name of the seafloor transponder positions (`fn1="pxp-ini.ini"` by default)
* `fn2`: Travel-time residual file name (`fn2="residual_sdls.out"` by default)
* `fno`: Output figure name (`fno="gradmap.pdf"` by default)
* `show`: if `show=true`, a figure is shown on REPL and is not saved as a file (`show=false` by default)
* `plot_size`: Figure size (`plot_size=(600,500)` by default)
* `lmargin`: Plot margin for the left edge (`lmargin=2.5` by default)
* `rmargin`: Plot margin for the right edge (`rmargin=1.0` by default)
* `tmargin`: Plot margin for the top edge (`tmargin=1.0` by default)
* `bmargin`: Plot margin for the bottom edge (`bmargin=1.0` by default)
* `ms1`: Plotted marker size for transponders (`ms1=10` by default)
* `ms2`: Plotted marker size for track (`ms2=3` by default)
* `gfs`: Label fontsize (`gfs=12` by default)
* `shape`: Plotted marker shape (`shape=:cross` by default)
* `cb`: Color scale (`cb=:lightrainbow` by default)

# Example
    plot_gradmap_gradv((-2500,2500),(-2500,2500),autoscale=false,fno="gradmap.png")

"""
function plot_gradmap_gradv(xrange=(-3000,3000),yrange=(-3000,3000);autoscale=true::Bool,fn1="pxp-ini.inp"::String,fn2="residual_sdls.out"::String,fno="gradmap.pdf"::String,plot_size=(600,500),lmargin=2.5,tmargin=1.0,bmargin=1.0,rmargin=2.5,show=false::Bool,ms1=10::Int64,ms2=3::Int64,gfs=12::Int64,shape=:cross,cb=:lightrainbow)
  a = DelimitedFiles.readdlm(fn1)
  numk = size(a)[1]
  px = a[1:numk,1]
  py = a[1:numk,2]
  a = DelimitedFiles.readdlm(fn2)
  x = a[:,3]
  y = a[:,4]
  d = (a[:,6] -a[:,8] - a[:,11])*1000
  dmax = maximum([abs(minimum(d)),abs(maximum(d))])
  dmin = -1*dmax
  p = scatter(px,py,markershape=:utriangle,markersize=ms1,c=RGB(0.8,0.8,0.8))
  if autoscale == true
    scatter!(p,x,y,marker_z=d,markershape=shape,c=cb,clims=(dmin,dmax),xlabel="Easting [m]",ylabel="Northing [m]",aspect_ratio = 1,legend=false,colorbar=true, framestyle=:box,colorbar_title="\n TT residuals [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin),markersize=ms2,guidefontsize=gfs)
  else
    scatter!(p,x,y,marker_z=d,markershape=shape,c=cb,clims=(dmin,dmax),xlabel="Easting [m]",ylabel="Northing [m]",xlims=xrange,ylims=yrange,aspect_ratio = 1,legend=false,colorbar=true,colorbar_title="\n TT residuals [msec]",framestyle=:box,size=plot_size,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin),markersize=ms2,guidefontsize=gfs)
  end
  if show == false
    savefig(p,fno)
  else
    gui(p)
  end
end

