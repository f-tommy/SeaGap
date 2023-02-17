#using Plots
#using Random
#using DelimitedFiles
#plot_mcmcparam()

export plot_mcmcparam
"""
    plot_mcmcparam(NA; fno,fn,plot_size,lmargin,rmargin,tmargin,bmargin,show,nshuffle,ms,lfs,tfs,bmargin0)

Make a figure of the unknown parameters change through the MCMC iteration using the results of `pos_array_mcmcgrad()` and `static_array_mcmcgradc()`.

* `NA`: Sampling interval of the MCMC prcessing (`NA=5` by default), which must be same with `NA` in `pos_array_mcmcpvg`
* `fn`: the input file name (`fn="sample.out"` by default)
* `fno`: Output figure name (`fno="mcmc_param.pdf"` by default)
* `show`: if `show=true`, a figure is shown on REPL and is not saved as a file (`show=false` by default)
* `nshuffle`: number of plots for each parameter (if all samples are plotted, the figure is crowded; thus, `nshuffle` of samples are randomly picked; if `nshuffle=0`, all samples are plotted; `nshuffle=10000` by default)
* `plot_size`: Figure size (`plot_size=(500,800)` by default)
* `lmargin`: Plot margin for the left edge (`lmargin=6.0` by default)
* `rmargin`: Plot margin for the right edge (`rmargin=1.0` by default)
* `tmargin`: Plot margin for the top edge (`tmargin=1.0` by default)
* `bmargin`: Plot margin for the bottom edge (`bmargin=1.0` by default)
* `bmargin0`: Plot margin for the bottom edges of upper two panels (`bmargin0=-2.0` by default)
* `ms`: Plotted marker size (`ms=2` by default)
* `lfs`: Label fontsize (`lfs=5` by default)
* `tfs`: Tick fontsize (`tfs=5` by default)

# Example
    plot_mcmcparam(5,fno="mcmc_param.pdf")
"""
function plot_mcmcparam(NA=5::Int64; fno="mcmc_param.pdf"::String,fn="sample.out"::String,plot_size=(500,800),lmargin=6,rmargin=1.0, tmargin=1.0, bmargin=1.0,show=false::Bool,nshuffle=10000::Int64,ms=2::Int64,lfs=5::Int64,tfs=5::Int64,bmargin0=-2)
  println(stderr," --- Read $fn")
  dat, list = DelimitedFiles.readdlm(fn,header=true)
  num = size(dat)[1]
  if num < 1
    error(" plot_mcmcparam: No data")
  end
  if nshuffle > 0 && nshuffle < num
    println(stderr," --- Resampling")
    nl = shuffle(Vector(1:num))[1:nshuffle]
    dx = zeros(nshuffle); dy = zeros(nshuffle); dz = zeros(nshuffle)
    gx = zeros(nshuffle); gy = zeros(nshuffle); gd = zeros(nshuffle)
    for n in 1:nshuffle
      dx[n] = dat[nl[n],1]
      dy[n] = dat[nl[n],2]
      dz[n] = dat[nl[n],3]
      gx[n] = dat[nl[n],4]*1000
      gy[n] = dat[nl[n],5]*1000
      gd[n] = dat[nl[n],6]
    end
  else
    dx = dat[:,1]
    dy = dat[:,2]
    dz = dat[:,3]
    gx = dat[:,4]*1000
    gy = dat[:,5]*1000
    gd = dat[:,6]
  end
  println(stderr," --- Plot")
  p0 = scatter(nl*NA,dx,ylabel="Displacement (EW) [m]",legend=:none,markershape=:cross,xtick=:none,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin0),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin),framestyle=:box,ms=ms)
  p1 = scatter(nl*NA,dy,ylabel="Displacement (NS) [m]",legend=:none,markershape=:cross,xtick=:none,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin0),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, rmargin),framestyle=:box,ms=ms)
  p2 = scatter(nl*NA,dz,ylabel="Displacement (UD) [m]",legend=:none,markershape=:cross,xtick=:none,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin0),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, rmargin),framestyle=:box,ms=ms)
  p3 = scatter(nl*NA,gx,ylabel="Shallow Gradient (EW) [msec/km]",legend=:none,markershape=:cross,xtick=:none,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin0),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, rmargin),framestyle=:box,ms=ms)
  p4 = scatter(nl*NA,gy,ylabel="Shallow Gradient (NS) [msec/km]",legend=:none,markershape=:cross,xtick=:none,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin0),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, rmargin),framestyle=:box,ms=ms)
  p5 = scatter(nl*NA,gd,ylabel="Gradient depth [km]",xlabel="MCMC Iterations after burn-in period",legend=:none,markershape=:cross,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, rmargin),framestyle=:box,ms=ms)
  plt = plot(p0,p1,p2,p3,p4,p5,layout=(6,1), size=plot_size,labelfontsize=lfs,tickfontsize=tfs)
  println(stderr," --- Save")
  if show == false
    savefig(plt,fno)
  else
    gui(plt)
  end
end
