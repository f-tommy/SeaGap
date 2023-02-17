#using Plots
#using DelimitedFiles
#using Statistics
#using Printf

export plot_denoise_kinematic, denoise_kinematic
"""
    plot_denoise_kinematic(;xrange,yrange,autoscale,fn1,fn2,fno,plot_size,lmargin1,lmargin2,rmargin,tmargin,bmargin,show,ms1,m2)

Make a figure of `fno` from the `denoise_kinematic()` output file.

* `xrange`: EW range for the array positions [m] (`xrange=(-1,1)` by default)
* `yrange`: NS range for the array positions [m] (`yrange=(-1,1)` by default)
* `autoscale`: if `autoscale=true`, the EW and NS ranges are automatically determined; if `autoscale=false`, the EW and NS ranges are fixed to `xrange` and `yrange` (`autoscale=true` by default)
* `plot_size`: Figure size (`plot_size=(900,600)` by default)
* `fn1`: Input file name (`fn="kinematic_array.out"` by defualt)
* `fn2`: Input file name after denoise (`fn="tmp"` by defualt)
* `fno`: Output file name (`fno="denoise_kinematic.pdf"` by defualt)
* `lmargin1`: Plot margin for the left edge for map (`lmargin=4.0` by default)
* `lmargin2`: Plot margin for the left edge for time-series (`lmargin=3.0` by default)
* `rmargin`: Plot margin for the right edge (`rmargin=1.0` by default)
* `tmargin`: Plot margin for the top edge (`tmargin=1.0` by default)
* `bmargin`: Plot margin for the bottom edge (`bmargin=3.0` by default)
* `show`: if `show=true`, a figure is temporally shown; if false, the figure is save as `fno` (`show=false` by default)
* `ms1`: Plotted marker size for map (`ms1=5` by default)
* `ms2`: Plotted marker size for time-series (`ms2=4` by default)
"""
function plot_denoise_kinematic(;xrange=(-1,1),yrange=(-1,1), autoscale=true::Bool,fno="denoise_kinematic.pdf"::String,fn1="kinematic_array.out"::String,fn2="tmp",plot_size=(900,600),lmargin1=4.0,lmargin2=4.0, rmargin=1.0, tmargin=1.0, bmargin=1.0, show=false::Bool, ms1=5::Int64,ms2=4::Int64)
  n, m, dat01 = read_matrix(fn1)
  dat1 = unixsort(dat01,1)
  n, m, dat02 = read_matrix(fn2)
  dat2 = unixsort(dat02,1)
  t01 = dat1[1,1]; t02 = dat2[1,1]
  t1 = (dat1[:,1] .- t01)./(60*60)
  t2 = (dat2[:,1] .- t02)./(60*60)
  if autoscale == true
    p1 = scatter(dat1[:,3],dat1[:,4],aspect_ratio=1,markershape=:cross,label=:none,mc=:red,framestyle=:box,markersize=ms1,xlabel="EW position [m]",ylabel="NS position [m]",left_margin=Plots.Measures.Length(:mm,lmargin1),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, 0.0))
    p2 = scatter(t1,dat1[:,3],markershape=:cross,label=:none,mc=:red,markersize=ms2,framestyle=:box,ylabel="EW position [m]",left_margin=Plots.Measures.Length(:mm, lmargin2),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin))
    p3 = scatter(t1,dat1[:,4],markershape=:cross,label=:none,mc=:red,framestyle=:box,markersize=ms2,xlabel="Time [hours]",ylabel="NS position [m]",left_margin=Plots.Measures.Length(:mm, lmargin2),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin))
  else
    p1 = scatter(dat1[:,3],dat1[:,4],framestyle=:box,aspect_ratio=1,xlim=xrange,ylim=yrange,markershape=:cross,label=:none,mc=:red,markersize=ms1,xlabel="EW position [m]",ylabel="NS position [m]",left_margin=Plots.Measures.Length(:mm, lmargin1),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, 0.0))
    p2 = scatter(t1,dat1[:,3],ylim=xrange,framestyle=:box,markershape=:cross,label=:none,mc=:red,markersize=ms2,ylabel="EW position [m]",left_margin=Plots.Measures.Length(:mm, lmargin2),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin))
    p3 = scatter(t1,dat1[:,4],ylim=yrange,markershape=:cross,framestyle=:box,label=:none,mc=:red,markersize=ms2,xlabel="Time [hours]",ylabel="NS position [m]",left_margin=Plots.Measures.Length(:mm, lmargin2),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin))
  end
  scatter!(p1,dat2[:,3],dat2[:,4],markershape=:cross,label=:none,mc=:blue,markersize=ms1)
  scatter!(p2,t2,dat2[:,3],markershape=:cross,label=:none,mc=:blue,markersize=ms2)
  scatter!(p3,t2,dat2[:,4],markershape=:cross,label=:none,mc=:blue,markersize=ms2)
  plts = plot(plot(p1),plot(p2,p3,layout=(2,1)),size=plot_size,layout=(1,2))
  if show == false
    savefig(plts,fno)
  else
    gui(plts)
  end
end

"""
    denoise_kinematic(;xrange,yrange,method,autoscale,n,sigma1,sigma2,type,save,prompt,fn,fn0,plot_size,lmargin1,lmargin2,rmargin,tmargin,bmargin,show,fno1,fno2)

Eliminate outliers from the results of `kinematic_array()`.
Calculate travel-time residual, estimate smoothed travel-time residuals by `n` running `method` filter, exclude outliers beyond `sigma` Std, and plot them by `plot_denoise()`. The denoised observational file is rewritten in `fn4`.

* `method`: Method of running filter ("mean" or "method"; `method="median"` by default)
* `n`: Window size for the running filter
* `sigma1`: Outlier threshold for "spatial"
* `sigma2`: Outlier threshold for "temporal"
* `type`: "spatial","temporal", or "both" filters are conducted (`type="both"` by default) 
* `save`: if `save=true`, the input data file `fn` is renamed and saved as `fn0` (`save=true` by default)
* `prompt`: if `prompt=true`, confirmation message is shown; if false, the input file is forcely saved (`prompt=true` by default)  

* `fn`: Input data file (`fn="kinematic_array.out"` by default)
* `fn0`: if `save=true`, `fn` is saved (`fn0="kinematic_array.out_org"` by default)
* `xrange`: EW range for the array positions [m] (`xrange=(-1,1)` by default)
* `yrange`: NS range for the array positions [m] (`yrange=(-1,1)` by default)
* `autoscale`: if `autoscale=true`, the vertical ranges are automatically determined; if `autoscale=false`, the EW and NS ranges are fixed to `xrange` and `yrange` (`autoscale=true` by default)
* `plot_size`: Figure size (`plot_size=(900,600)` by default)
* `lmargin1`: Plot margin for the left edge of map (`lmargin=4.0` by default)
* `lmargin2`: Plot margin for the left edge of time-series (`lmargin=3.0` by default)
* `rmargin`: Plot margin for the right edge (`rmargin=1.0` by default)
* `tmargin`: Plot margin for the top edge (`tmargin=1.0` by default)
* `bmargin`: Plot margin for the bottom edge (`bmargin=3.0` by default)
* `show`: if `show=true`, a figure is temporally shown; if false, the figure is save as `fno` (`show=true` by default)
* `ms1`: Plotted marker size for map (`ms1=5` by default)
* `ms2`: Plotted marker size for time-series (`ms2=4` by default)

* `fno1`: Output figure name
* `fno2`: Output text file for the elimited data list

# Example
    denoise_kinematic(n=7,sigma1=4.0,sigma2=5.0,method="median")
"""
function denoise_kinematic(;xrange=(-1,1),yrange=(-1,1),method="median"::String,autoscale=true::Bool,n=15,sigma1=4.0,sigma2=4.0,type="both",save=true::Bool,prompt=true::Bool,fn="kinematic_array.out"::String, fn0="kinematic_array.out_org"::String, plot_size=(1000,500),lmargin1=4.0,lmargin2=3.0,rmargin=1.0, tmargin=1.0, bmargin=3.0, show=true::Bool, fno1="kinematic_denoise.pdf"::String,ms1=5::Int64,ms2=4::Int64,fno2="eliminated_list.out")
  # Read file
  N, M, dat0 = read_matrix(fn)
  dat = unixsort(dat0,1)
  x0 = dat[:,3]; y0 = dat[:,4]
  if save ==true 
    cp("$fn","$fn0",force=true)
    println(stderr,"Copy $fn -> $fn0")
  end
  dlist = []
  # Spatial filter
  if type != "temporal"
    println(stderr,"  Perform spatial filtering")
    if method == "median"
      ax = median(x0); ay = median(y0)
    else
      ax = mean(x0); ay = mean(y0)
    end
    xv = x0 .- ax; yv = y0 .- ay
    xr = std(xv); yr = std(yv)
    x1 = ax - sigma1*xr; x2 = ax + sigma1*xr
    y1 = ay - sigma1*yr; y2 = ay + sigma1*yr
    dat1 = dat[dat[:,3].>x1,:]
    dlist = vcat(dlist,dat[dat[:,3].<=x1,9])
    dat2 = dat1[dat1[:,3].<x2,:]
    dlist = vcat(dlist,dat1[dat1[:,3].>=x2,9])
    dat3 = dat2[dat2[:,4].<y2,:]
    dlist = vcat(dlist,dat2[dat2[:,4].>=y2,9])
    dat = dat3[dat3[:,4].>y1,:]
    dlist = vcat(dlist,dat3[dat3[:,4].<=y1,9])
  end
  x0 = dat[:,3]; y0 = dat[:,4]
  # Running filter
  if type != "spatial"
    println(stderr,"  Perform temporal filtering")
    if method == "median"
      xm = runmed(x0,n)
      ym = runmed(y0,n)
    else
      xm = runave(x0,n)
      ym = runave(y0,n)
    end
    dx = x0 - xm
    dy = y0 - ym
    sx = sigma2*std(dx)
    sy = sigma2*std(dy)
    dat0 = hcat(dat,abs.(dx),abs.(dy))
    dat1 = dat0[dat0[:,10].<sx,:]
    dlist = vcat(dlist,dat0[dat0[:,10].>=sx,9])
    dat = dat1[dat1[:,11].<sy,:]
    dlist = vcat(dlist,dat1[dat1[:,11].>=sy,9])
  end
  # Make a new observational data file
  println(stderr,"Make denoised position data file: $fn")
  open("tmp","w") do out
    num = size(dat)[1]
    for i in 1:num
      @printf(out,"%10.4f %i %4.6f %4.6f %.6e %s %s %s %i\n",dat[i,1],dat[i,2],dat[i,3],dat[i,4],dat[i,5],dat[i,6],dat[i,7],dat[i,8],dat[i,9])
    end
  end
  # Plot
  println(stderr,"Plot denoised time-series")
  plot_denoise_kinematic(fno=fno1,fn1=fn,fn2="tmp",xrange=xrange,yrange=yrange,plot_size=plot_size,lmargin1=lmargin1,lmargin2=lmargin2,rmargin=rmargin, tmargin=tmargin, bmargin=bmargin, show=show,ms1=ms2,ms2=ms2) 
  println(stderr,dlist)
  if prompt == true
    q = Base.prompt("Do you accept the denoise processing? (yes/no)")
    if q == "yes"
      mv("tmp",fn,force=true)
      ne = size(dlist)[1]
      open(fno2,"a") do out2
        for n in 1:ne
          @printf(out2,"%i\n",dlist[n])
        end
      end
    else
      rm("tmp")
    end
  else
    mv("tmp",fn,force=true)
    ne = size(dlist)[1]
    open(fno2,"a") do out2
      for n in 1:ne
        @printf(out2,"%i\n",dlist[n])
      end
    end
  end
end 
