#using Plots
#using DelimitedFiles

export plot_gradmap
function plot_gradmap(xrange=(-3000,3000),yrange=(-3000,3000);autoscale=true::Bool,fn1="pxp-ini.xyh"::String,fn2="residual_grad.out"::String,fno="gradmap.pdf"::String,plot_size=(600,500),lmargin=2.5,tmargin=1.0,bmargin=1.0,rmargin=1.0,show=false::Bool,ms1=10::Int64,ms2=3::Int64,gfs=12::Int64,shape=:circle,cb=:lightrainbow)
  a = DelimitedFiles.readdlm(fn1)
  numk = size(a)[1]
  px = a[1:numk,1]
  py = a[1:numk,2]
  a = DelimitedFiles.readdlm(fn2)
  x = a[:,3]
  y = a[:,4]
  d = (a[:,6] -a[:,8] - a[:,10])*1000
  dmax = maximum([abs(minimum(d)),abs(maximum(d))])
  dmin = -1*dmax
  p = scatter(px,py,markershape=:utriangle,markersize=ms1,c=RGB(0.8,0.8,0.8))
  if autoscale == true
    scatter!(p,x,y,marker_z=d,markershape=shape,c=cb,clims=(dmin,dmax),xlabel="Easting [m]",ylabel="Northing [m]",aspect_ratio = 1,legend=false,colorbar=true, framestyle=:box,colorbar_title="TT residuals [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin),markersize=ms2,guidefontsize=gfs)
  else
    scatter!(p,x,y,marker_z=d,markershape=shape,c=cb,clims=(dmin,dmax),xlabel="Easting [m]",ylabel="Northing [m]",xlims=xrange,ylims=yrange,aspect_ratio = 1,legend=false,colorbar=true,colorbar_title="TT residuals [msec]",framestyle=:box,size=plot_size,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin),markersize=ms2,guidefontsize=gfs)
  end
  if show == false
    savefig(p,fno)
  else
    gui(p)
  end
end

