#using Plots
#using DelimitedFiles

export plot_time_array_each, plot_map_array_each
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
    p2 = scatter(t,x,xlabel="Time [hour]",ylabel="Northing [m]",framestyle=:box,legend = :none,markershape=:cross,bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, 0),markersize=ms)
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


