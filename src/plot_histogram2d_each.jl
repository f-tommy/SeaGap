#using Dates
#using DelimitedFiles
#using Plots

export plot_histogram2d_each
"""
    plot_histogram2d_each(param1,param2; fn,show,fno,plot_size,lmargin,tmargin,bmargin,rmargin,nbins)

Make a figure of a heatmap for certain two parameters (`param1` and `param2`) from the sampling results of `static_array_mcmcgrad()` or `static_array_mcmcgradc()`.

* `param1`: Parameter name 1
* `param2`: Parameter name 2
* `fn`: Input file name (`fn="sample.out"` by default)
* `fno`: Output figure name (`fno="histogram2d_each.pdf"` by default)
* `plot_size`: Figure size (`plot_size=(650,600)` by default)
* `lmargin`: Plot margin for the left edge (`lmargin=1.5` by default)
* `rmargin`: Plot margin for the right edge (`rmargin=1.5` by default)
* `tmargin`: Plot margin for the top edge (`tmargin=0.5` by default)
* `bmargin`: Plot margin for the bottom edge (`bmargin=0.5` by default)
* `show`: if `show=true`, a figure is temporally shown; if false, the figure is saved as `fno` (`show=false` by default)
* nbins: Number of histogram's intervals (`nbins=50` by default)

# Example
    plot_histogram2d_each("UD_disp.","S-NTD_10",fno="histogram2d_UD_SNTD-10.png")
"""
function plot_histogram2d_each(param1="EW_disp."::String,param2="NS_disp."::String; fn="sample.out"::String,show=false::Bool,fno="histogram2d_each.pdf"::String,plot_size=(650,600),lmargin=1.5,tmargin=0.5,bmargin=0.5,rmargin=1.5,nbins=50::Int64)
  println(stderr," === Drawing 2d-histogram for static_array_mcmcgrad samples ===")
  time1 = now()
  # --- Read data
  println(stderr," --- Read files")
  dat, list = DelimitedFiles.readdlm(fn, header=true)
  num = length(list)
  num1 = 0; num2 = 0
  for n in 1:num
    if list[n] == param1
      num1 = n
    end
    if list[n] == param2
      num2 = n
    end
  end
  if num1 == 0
    println(stderr,list)
    error("$param1 is not properly given")
  end
  if num2 == 0
    println(stderr,list)
    error("$param2 is not properly given")
  end
  # --- Plot
  println(stderr,"--- Drawing histogram $num1 $num2")
  plt = histogram2d(dat[:,num1],dat[:,num2],nbins=nbins,xlabel=param1,ylabel=param2,c=:dense,framestyle=:box,size=plot_size,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin))
  if show == false
    savefig(plt,fno)
  else
    gui(plt)
  end
  # --- Close process
  time2 = now()
  println(stderr," Start time:",time1)
  println(stderr," Finish time:",time2)
end
