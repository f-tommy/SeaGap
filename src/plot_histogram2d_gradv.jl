#using Dates
#using Statistics
#using DelimitedFiles
#using Random
#using Plots
#plot_histogram2d()

export plot_histogram2d_gradv
"""
    plot_histogram2d_gradv(;fn,show,fno,plot_size,lmargin,tmargin,bmargin,rmargin,nshuffle,lfs,tfs,hm,vm,nbins,pscale)

Make a figure showing histograms, heatmaps, and scatter maps for major seven parameters (array displacements, shallow gradients, gradient depth).

* `fn`: Input file (`fn="sample.out"` by default)
* `fno`: Output figure name
* `show`: if `show=true`, a figure is shown on REPL and is not saved as a file (`show=false` by default)
* `nshuffle`: number of plots for each parameter (if all samples are plotted, the figure is crowded; thus, `nshuffle` of samples are randomly picked; if `nshuffle=0`, all samples are plotted; `nshuffle=10000` by default)
* `nbins`: Number of histogram's intervals (`nbins=50` by default)
* `plot_size`: Figure size (`plot_size=(700,700)` by default)
* `lmargin`: Plot margin for the left edge (`lmargin=1.5` by default)
* `rmargin`: Plot margin for the right edge (`rmargin=1.5` by default)
* `tmargin`: Plot margin for the top edge (`tmargin=1.0` by default)
* `bmargin`: Plot margin for the bottom edge (`bmargin=1.0` by default)
* `lfs`: Label fontsize (`lfs=5` by default)
* `tfs`: Tick fontsize (`tfs=4` by default)
* `pscale`: Length scaling factor for horizontal axis (`pscale=3` by default)
* `hm`: Horizontal interval of panels (`hm=-1.0` by default) 
* `vm`: Vertical interval of panels (`vm=-1.0` by default) 

# Example
    plot_histogram2d_gradv(fno="histogram2d.pdf")
"""
function plot_histogram2d_gradv(;fn="sample.out"::String,show=false::Bool,fno="histogram2d.pdf"::String,plot_size=(750,750),lmargin=1.5,tmargin=1.0,bmargin=1.0,rmargin=1.5,nshuffle=10000::Int64,lfs=5::Int64,tfs=4::Int64,hm=-1.0,vm=-1.0,nbins=50::Int64,pscale=3)
  println(stderr," === Drawing 2d-histograms for static_array_mcmcgradv samples ===")
  time1 = now()
  # --- Read data
  println(stderr," --- Read files")
  dat0, list0 = DelimitedFiles.readdlm(fn, header=true)
  dat = hcat(dat0[:,1:3]*100,dat0[:,4:5]*1000,dat0[:,6:7])
  list = ["Disp. (EW) [cm]","Disp. (NS) [cm]","Disp. (UD) [cm]","S-Grad. (EW)\n [msec/km]","S-Grad. (NS)\n [msec/km]","G-Depth (EW)\n [km]","G-Depth (NS)\n [km]"]
  num = size(dat)[1]
  if num < 1
    error(" plot_histogram2d: No data")
  end
  # --- Resampling
  if nshuffle > 0
    println(stderr," --- Resampling for upper-right side scatter plots")
    nl = shuffle(Vector(1:num))[1:nshuffle]
    datr = zeros(nshuffle,7)
    for n in 1:nshuffle
      for m = 1:7
        datr[n,m] = dat[nl[n],m]
      end
    end
  else
    datr = copy(dat)
  end
  # --- Statistical setting
  dmin = minimum(dat,dims=1)
  dmax = maximum(dat,dims=1)
  dt = (dmax - dmin) / pscale
  lmin = dmin - dt
  lmax = dmax + dt
  # --- Plot
  println(stderr," --- Drawing histograms")
  plt = []
  for j in 1:7
    for i in 1:7
      println(stderr,"   Plot: $i $j")
      if i == j
        if i == 7
          push!(plt, histogram(dat[:,i],bins=range(dmin[i],dmax[i],length=50),normed = true,xlims=(lmin[i],lmax[i]),xlabel=list[i],yticks=false,lc=:white,fc=:skyblue,legend=:none,framestyle=:box,lw=0.1,grid=false))
        else
          push!(plt, histogram(dat[:,i],bins=range(dmin[i],dmax[i],length=50),normed=true,xlims=(lmin[i],lmax[i]),lc=:white,fc=:skyblue,xticks=false,yticks=false,legend=:none,framestyle=:box,grid=false,lw=0.1))
        end
      elseif i < j
        if i == 1
          if  j == 7
            push!(plt, histogram2d(dat[:,i],dat[:,j],nbins=nbins,xlims=(lmin[i],lmax[i]),ylims=(lmin[j],lmax[j]),xlabel=list[i],ylabel=list[j],c=:dense,framestyle=:box,colorbar=:none,grid=false))
          else
            push!(plt, histogram2d(dat[:,i],dat[:,j],nbins=nbins,xlims=(lmin[i],lmax[i]),ylims=(lmin[j],lmax[j]),ylabel=list[j],xticks=false,c=:dense,framestyle=:box,colorbar=:none,grid=false))
          end
        else
          if  j == 7
            push!(plt, histogram2d(dat[:,i],dat[:,j],nbins=nbins,xlims=(lmin[i],lmax[i]),ylims=(lmin[j],lmax[j]),xlabel=list[i],yticks=false,c=:dense,framestyle=:box,colorbar=:none,grid=false))
          else
            push!(plt, histogram2d(dat[:,i],dat[:,j],nbins=nbins,xlims=(lmin[i],lmax[i]),ylims=(lmin[j],lmax[j]),xticks=false,yticks=false,grid=false,c=:dense,framestyle=:box,colorbar=:none))
          end
        end
      else
        if i == 7
          push!(plt, scatter(datr[:,i],datr[:,j],xlims=(lmin[i],lmax[i]),ylims=(lmin[j],lmax[j]),xticks=false,grid=false,ylabel=list[j],framestyle=:box,mc=:skyblue,ms=1,marker=:circle,msw=0,ymirror=true,legend=:none))
        else
          push!(plt, scatter(datr[:,i],datr[:,j],xlims=(lmin[i],lmax[i]),ylims=(lmin[j],lmax[j]),xticks=false,yticks=false,grid=false,framestyle=:box,mc=:skyblue,ms=1,marker=:circle,msw=0,legend=:none))
        end
      end
    end
  end
  println(stderr,"   Layouting")
  plts = plot(plt...,layout=(7,7),size=plot_size,labelfontsize=lfs,tickfontsize=tfs,left_margin=Plots.Measures.Length(:mm,hm),bottom_margin=Plots.Measures.Length(:mm,vm),top_margin=Plots.Measures.Length(:mm,vm),right_margin=Plots.Measures.Length(:mm,hm))
  if show == false
    savefig(plts,fno)
  else
    gui(plts)
  end
  # --- Close process
  time2 = now()
  println(stderr," Start time:",time1)
  println(stderr," Finish time:",time2)
end
