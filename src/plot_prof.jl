#using Plots
#using DelimitedFiles

export plot_prof
"""
    plot_prof(;fno,fn,plot_size,lmargin,tmargin,bmargin,rmargin,show)

Make a figure of a sound speed profile.

* `fn`: Input file of a sound speed profile (`fn="ss_prof.zv"` in default)
* `fno`: Output figure name (`fno="ss_prof.pdf"` in default)
* `plot_size`: Figure size (`plot_size=(450,600)` in default)
* `show`: if `show=true`, a figure is shown on REPL and is not saved as a file (`show=false` in default)
* `lmargin`: Plot margin for the left edge (`lmargin=2.5` in default)
* `rmargin`: Plot margin for the right edge (`rmargin=1.0` in default)
* `tmargin`: Plot margin for the top edge (`tmargin=1.0` in default)
* `bmargin`: Plot margin for the bottom edge (`bmargin=1.0` in default)

# Example
    plot_prof(fno="ss_prof.png",fn="ss_prof.zv",show=false)
"""
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
