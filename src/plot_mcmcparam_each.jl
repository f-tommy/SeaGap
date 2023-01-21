#using Plots
#using Random
#using DelimitedFiles
# Usage: plot_mcmcparam_each("S-NTD_10")

export plot_mcmcparam_each
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
