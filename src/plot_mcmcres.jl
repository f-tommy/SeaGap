#using Plots
#using Random
#using DelimitedFiles
#plot_mcmcres()

export plot_mcmcres
"""
    plot_mcmcres(;fno,fn,plot_size,lmargin,rmargin,tmargin,bmargin,show,nshuffle,ms,lfs,tfs,bmargin0)

Make a figure of the posterior PDF and the travel-time RMS changes through the MCMC iteration using the results of `static_array_mcmcgrad()` and `static_array_mcmcgradc()`.

* `fn`: the input file name (`fn="mcmc.out"` by default)
* `fno`: Output figure name (`fno="mcmc_res.pdf"` by default)
* `show`: if `show=true`, a figure is shown on REPL and is not saved as a file (`show=false` by default)
* `nshuffle`: number of plots (if all samples are plotted, the figure is crowded; thus, `nshuffle` of samples are randomly picked; if `nshuffle=0`, all samples are plotted; `nshuffle=50000` by default)
* `plot_size`: Figure size (`plot_size=(600,800)` by default)
* `lmargin`: Plot margin for the left edge (`lmargin=3.0` by default)
* `rmargin`: Plot margin for the right edge (`rmargin=1.0` by default)
* `tmargin`: Plot margin for the top edge (`tmargin=1.0` by default)
* `bmargin`: Plot margin for the bottom edge (`bmargin=1.0` by default)
* `bmargin0`: Plot margin for the bottom edges of upper two panels (`bmargin0=-1.0` by default)
* `ms`: Plotted marker size (`ms=3` by default)
* `lfs`: Label fontsize (`lfs=10` by default)
* `tfs`: Tick fontsize (`tfs=6` by default)

# Example
    plot_mcmcres(fno="mcmc_res.pdf")
"""
function plot_mcmcres(;fno="mcmc_res.pdf"::String,fn="mcmc.out"::String,plot_size=(600,400),lmargin=3.0,rmargin=1.0, tmargin=1.0, bmargin=1.0,show=false,nshuffle=50000::Int64,ms=3::Int64,lfs=10::Int64,tfs=6::Int64,bmargin0=-1)
  println(stderr," --- Read $fn")
  dat = DelimitedFiles.readdlm(fn)
  num = size(dat)[1]
  if num < 1
    error(" plot_mcmcres: No data")
  end
  if nshuffle > 0 && nshuffle < num
    println(stderr," --- Resampling")
    nl = shuffle(Vector(1:num))[1:nshuffle]
    t = zeros(nshuffle); h = zeros(nshuffle); r = zeros(nshuffle)
    for n in 1:nshuffle
      t[n] = dat[nl[n],1]
      h[n] = dat[nl[n],6]
      r[n] = dat[nl[n],7]*1000
    end
  else
    t = dat[:,1]
    h = dat[:,6]
    r = dat[:,7]*1000
  end
  println(stderr," --- Plot")
  p0 = scatter(t,h,ylabel="Posterior PDF",legend=:none,markershape=:cross,xtick=:none,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin0),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin),framestyle=:box,ms=ms)
  p1 = scatter(t,r,xlabel="MCMC Iteration",ylabel="TT RMS [msec]",legend=:none,markershape=:cross,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, rmargin),framestyle=:box,ms=ms)
  plt = plot(p0,p1,layout=(2,1), size=plot_size, labelfontsize=lfs, tickfontsize=tfs)
  println(stderr," --- Save")
  if show == false
    savefig(plt,fno)
  else
    gui(plt)
  end
end
