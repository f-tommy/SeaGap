#using Dates
#using DelimitedFiles
#using Plots

export plot_histogram2d_each
function plot_histogram2d_each(param1="EW_disp."::String,param2="NS_disp."::String; fn="sample.out"::String,show=false::Bool,fno="histogram2d_each.pdf"::String,plot_size=(650,600),lmargin=1.5,tmargin=0.5,bmargin=0.5,rmargin=1.5,nbins=50::Int64)
  println(stderr," === Drawing 2d-histogram for pos_array_mcmcpvg samples ===")
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
