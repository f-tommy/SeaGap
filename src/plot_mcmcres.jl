#using Plots
#using Random
#using DelimitedFiles
#plot_mcmcres()

export plot_mcmcres
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
