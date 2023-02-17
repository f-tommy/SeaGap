#using Plots
#using Random
#using DelimitedFiles
# Usage: plot_mcmcparam_each("S-NTD_10")

export plot_mcmcparam_each
"""
    plot_mcmcparam_each(param,NA; fno,fn,plot_size,lmargin,rmargin,tmargin,bmargin,show,nshuffle,ms,lfs,tfs)

Make a figure of an unknown parameter change through the MCMC iteration using the results of `static_array_mcmcgrad()` and `static_array_mcmcgradc()`.

* `param`: Parameter name
* `NA`: Sampling interval of the MCMC prcessing (`NA=5` by default), which must be same with `NA` in `pos_array_mcmcpvg`
* `fn`: the input file name (`fn="sample.out"` by default)
* `fno`: Output figure name (`fno="mcmc_param_each.pdf"` by default)
* `show`: if `show=true`, a figure is shown on REPL and is not saved as a file (`show=false` by default)
* `nshuffle`: number of plots for each parameter (if all samples are plotted, the figure is crowded; thus, `nshuffle` of samples are randomly picked; if `nshuffle=0`, all samples are plotted; `nshuffle=10000` by default)
* `plot_size`: Figure size (`plot_size=(500,200)` by default)
* `lmargin`: Plot margin for the left edge (`lmargin=2.0` by default)
* `rmargin`: Plot margin for the right edge (`rmargin=1.0` by default)
* `tmargin`: Plot margin for the top edge (`tmargin=1.0` by default)
* `bmargin`: Plot margin for the bottom edge (`bmargin=1.5` by default)
* `ms`: Plotted marker size (`ms=2` by default)
* `lfs`: Label fontsize (`lfs=7` by default)
* `tfs`: Tick fontsize (`tfs=5` by default)

# Example
    plot_mcmcparam_each("S-NTD_10",5,fno="mcmc_param_S-NTD_10.png")
"""
function plot_mcmcparam_each(param::String,NA=5::Int64; fno="mcmc_param_each.pdf"::String,fn="sample.out"::String,plot_size=(500,200),lmargin=2,rmargin=1.0, tmargin=1.0, bmargin=1.5,show=false,nshuffle=20000::Int64,ms=2::Int64,lfs=7::Int64,tfs=5::Int64)
  println(stderr," --- Read $fn")
  dat, list = DelimitedFiles.readdlm(fn,header=true)
  num, mm = size(dat)
  if num < 1
    error(" plot_mcmcparam_each: No data")
  end
  nn = 0
  for m in 1:mm
    if param == list[m]
      nn = m
    end
  end
  if nshuffle > 0 && nshuffle < num
    println(stderr," --- Resampling")
    nl = shuffle(Vector(1:num))[1:nshuffle]
    d = zeros(nshuffle)
    for n in 1:nshuffle
      d[n] = dat[nl[n],nn]
    end
  else
    d = dat[:,nn]
  end
  println(stderr," --- Plot")
  plt = scatter(nl*NA,d,ylabel=param,xlabel="MCMC Iterations after burn-in period",legend=:none,markershape=:cross,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin),framestyle=:box,ms=ms,size=plot_size,labelfontsize=lfs,tickfontsize=tfs)
  println(stderr," --- Save")
  if show == false
    savefig(plt,fno)
  else
    gui(plt)
  end
end
