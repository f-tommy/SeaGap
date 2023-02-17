#using Dates
#using Statistics
#using DelimitedFiles
#using Plots
# Usage: plot_cormap()
# Usage: plot_cormap(type="all",NPB=41,as=5,fno="cormap_all.pdf")

export plot_cormap

"""
    plot_cormap(;txt,all,fn,fno,fno0,as,pfs,plot_size,show,lmargin,tmargin,bmargin,rmargin)

Make a figure correlation coefficient map from the sampling results of `static_array_mcmcgrad()` or `static_array_mcmcgradc()`.

* `txt`: if `txt=true`, text data file of correlation coefficients is saved (`txt=false` by default)
* `all`: if `all=true`, correlation coefficients are calculated for all parameters; if `all=false`, correlation coefficients are calculated for major paramters
* `fn`: the input file name (`fn="sample.out"` by default)
* `fno0`: the output text file name (`fno0="correlation.out"` by default)
* `fno`: the output figure file name (`fno="cormap.pdf"` by default)
* `as`: Tick (parameter name) fontsize (`ts=10` by default)
* `pfs`: Annotation (correlation coefficients) fontsize (i`pfs=8` by default)
* `plot_size`: Figure size (`plot_size=(600,600)` by default)
* `lmargin`: Plot margin for the left edge (`lmargin=0.5` by default)
* `rmargin`: Plot margin for the right edge (`rmargin=4.0` by default)
* `tmargin`: Plot margin for the top edge (`tmargin=0.5` by default)
* `bmargin`: Plot margin for the bottom edge (`bmargin=0.5` by default)
* `show`: if `show=true`, a figure is temporally shown; if false, the figure is save as `fno` (`show=false` by default)

# Example
    plot_cormap(type="all",NPB=41,as=5,fno="cormap_all.pdf")

"""
function plot_cormap(;txt=false::Bool,all=false::Bool,fn="sample.out"::String,fno="cormap.pdf"::String,fno0="correlation.out"::String,as=10::Int64,pfs=8::Int64,plot_size=(600,600),show=false::Bool,lmargin=0.5,tmargin=0.5,bmargin=0.5,rmargin=4.0)
  println(stderr," === Correlation for static_array_mcmcgrad samples ===")
  time1 = now()
  # --- Read data
  println(stderr," --- Read files")
  dat0, list0 = DelimitedFiles.readdlm(fn, header=true)
  if all == false
    dat = dat0[1:end,1:6]
    list = list0[1:6]
  else
    dat = copy(dat0)
    list = list0[1:end]
  end
  # --- Correlation
  println(stderr," --- Correlation")
  cdat = cor(dat)
  # --- Plot
  n,m = size(cdat)
  plt = heatmap(cdat,aspect_ratio=1,xtickfontsize=as,ytickfontsize=as,framestyle=:box,c=cgrad(:bwr),clims=(-1,1),xticks=(1:m,list), xrot=90, yticks=(1:m,list), yflip=true,size=plot_size,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin))
  if all == false
    annotate!(plt,[(j, i, text(round(cdat[i,j],digits=2), pfs,"Computer Modern",:black)) for i in 1:n for j in 1:m])
  end
  # --- plot
  if show == false
    savefig(plt,fno)
  else
    gui(plt)
  end
  # --- Output
  if txt == true
    println(stderr," --- Output")
    open(fno0,"w") do out
      print(out,"Name: ")
      for i in 1:m
        print(out,list[i]," ")
      end
      println(out,"")
      Base.print_array(out,hcat(list,cdat))
      println(out,"")
    end
  end
  # --- Close process
  time2 = now()
  println(stderr," Start time:",time1)
  println(stderr," Finish time:",time2)
end
