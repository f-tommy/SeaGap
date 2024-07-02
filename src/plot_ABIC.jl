#using Plots
#using DelimitedFiles

export plot_ABIC
"""
    plot_ABIC(;type,fno,fn,lmargin,rmargin,tmargin,bmargin,show)

Make a figure plotting the ABIC values searched by `static_array_s_ABIC()`.

* `fno`: Output figure name (`fno="ABIC_search.pdf"` by default)
* `fn`: Input file name obtained by `static_array_s_ABIC()` (`fn="ABIC_search.out"` by default)
* `lmargin`: Margin of figure for left edge (`lmargin=3.5` by default)
* `rmargin`: Margin of figure for right edge (`rmargin=1.5` by default)
* `tmargin`: Margin of figure for top edge (`tmargin=1.5` by default)
* `bmargin`: Margin of figure for bottom edge (`bmargin=1.5` by default)
* `show`: if `show=true`, a figure is temporally shown; if false, the figure is saved as `fno` (`show=falsei` by default)

# Example
    plot_ABIC()
"""
function plot_ABIC(;fno="ABIC_search.pdf"::String,fn="ABIC_search.out"::String, lmargin=3.5, rmargin=1.5, tmargin=1.5, bmargin=1.5, show=false::Bool)
  dat0 = DelimitedFiles.readdlm(fn)
  num = size(dat0)[1]
  dat = unixsort(dat0,1)
  k = dat[1:num,1]
  a = dat[1:num,5]
  minv, id = findmin(a)
  maxv = findmax(a)[1]
  mink = k[id]
  dv = (maxv-minv)/10
  ulim = maxv + dv
  llim = minv - dv
  p = plot([mink],st=:vline,legend=:none,c=:red,linewidth=1,ylim=(llim,ulim))
  plot!(p,annotations=(mink,minv, ("$mink", 8, 0.0, :top)))
  plot!(p,k,a,lc=:blue,framestyle=:box,xlabel="Hyper-parameter for the S-NTD norm",ylabel="ABIC",legend = :none,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin))
  scatter!(p,k,a,mc=:blue,legend = :none,markershape=:cross)
  if show == false
    savefig(p,fno)
  else
    gui(p)
  end
end
