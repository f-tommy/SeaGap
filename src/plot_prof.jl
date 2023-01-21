#using Plots
#using DelimitedFiles

export plot_prof
function plot_prof(;fno="ss_prof.pdf"::String,fn="ss_prof.zv"::String,plot_size=(450,600),lmargin=2.5,tmargin=1.0,bmargin=1.0,rmargin=1.0,show=false::Bool)
  a = DelimitedFiles.readdlm(fn)
  numz = size(a)[1]
  z = a[1:numz,1]
  v = a[1:numz,2]
  p = plot(v,z,xlabel="Sound speed [m/s]",ylabel="Water depth [m]",yflip=true,legend = :none,framestyle=:box,size=plot_size,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin))
  if show == false
    savefig(p,fno)
  else
    gui(p)
  end
end
