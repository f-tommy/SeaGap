#using Plots
#using DelimitedFiles

export plot_AICBIC
function plot_AICBIC(;type="BIC"::String,fno="AICBIC_search.pdf"::String,fn="AICBIC_search.out"::String, lmargin=3.5, rmargin=1.5, tmargin=1.5, bmargin=1.5, show=false::Bool)
  dat0 = DelimitedFiles.readdlm(fn)
  num = size(dat0)[1]
  dat = unixsort(dat0,1)
  k = dat[1:num,1]
  a = dat[1:num,2]
  b = dat[1:num,3]
  if type == "AIC"
    minv, id = findmin(a)
    maxv = findmax(a)[1]
    mink = Int64(floor(k[id]))
    dv = (maxv-minv)/10
    ulim = maxv + dv
    llim = minv - dv
    p = plot([mink],st=:vline,legend=:none,c=:red,linewidth=1,ylim=(llim,ulim))
    plot!(p,annotations=(mink,minv, ("$mink", 8, 0.0, :top)))
    plot!(p,k,a,lc=:blue,framestyle=:box,xlabel="Number of B-spline bases",ylabel=type,legend = :none,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin))
    scatter!(p,k,a,mc=:blue,legend = :none,markershape=:cross)
  else
    minv, id = findmin(b)
    maxv = findmax(b)[1]
    mink = Int64(floor(k[id]))
    dv = (maxv-minv)/10
    ulim = maxv + dv
    llim = minv - dv
    p = plot([mink],st=:vline,legend=:none,c=:red,linewidth=1,ylim=(llim,ulim))
    plot!(p,annotations=(mink,minv, ("$mink", 8, 0.0, :top)))
    plot!(p,k,b,lc=:blue,framestyle=:box,xlabel="Number of B-spline bases",ylabel=type,legend = :none,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin))
    scatter!(p,k,b,mc=:blue,legend = :none,markershape=:cross)
  end
  if show == false
    savefig(p,fno)
  else
    gui(p)
  end
end
