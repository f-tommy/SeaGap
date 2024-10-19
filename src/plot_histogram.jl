#using Dates
#using Statistics
#using DelimitedFiles
#using Plots
#using PDFmerger
# Usage: plot_histogram(drawnls=true)
# Usage: plot_histogram(NPB,all=true,fno="histogram_all.pdf")
# Usage: plot_histogram(NPB,all=true,fno="histogram_allv.pdf",drawnls=true)

export plot_histogram
"""
    plot_histogram(NPB; all,drawnls,fn,fn0,fno,plot_size,lmargin,tmargin,bmargin,rmargin,nbins)

Make a figure file of histograms from the sampling results obtained by `static_array_mcmcgrad()` or `static_array_mcmcgradc()`.

* `NPB`: Number of 3d B-spline bases
* `fn0`: Inversion results by `pos_array_all()` ("solve.out" by default), which is used when you'd like to plot the `static_array()` results on the histogram
* `fn`: Input file ("sample.out" by default)
* `fno`: Output figure name (note that this file must be a PDF file: "histogram.pdf" by default)
* `all`: if `all=true`, histograms for all parameters are drawn; if `all=false`, histograms for major six parameters (array displacements, shallow gradients, gradient depth) (`all=false` by default)
* `drawnls`: if `drawnls=true`, a normal distribution estimated by `static_array()` is drawn in a histogram (`drawnls=false` by default); the normal distributions are shown for the array displacements and 3d B-spline NTDs
* `nbins`: Number of histogram's intervals (`nbins=50` by default)
* `plot_size`: Figure size (`plot_size=(650,500)` by default)
* `lmargin`: Plot margin for the left edge (`lmargin=1.5` by default)
* `rmargin`: Plot margin for the right edge (`rmargin=1.5` by default)
* `tmargin`: Plot margin for the top edge (`tmargin=1.0` by default)
* `bmargin`: Plot margin for the bottom edge (`bmargin=1.0` by default)

# Example
    plot_histogram(77,all=true,fno="histogram_all.pdf")

"""
function plot_histogram(NPB=100::Int64; all=false,drawnls=false,fn="sample.out"::String,fn0="solve.out"::String,fno="histogram.pdf"::String,plot_size=(650,500),lmargin=1.5,tmargin=1.0,bmargin=1.0,rmargin=1.5,nbins=50::Int64)
  println(stderr," === Drawing histograms for static_array_mcmcgrad samples ===")
  time1 = now()
  println(stderr," --- Output file check")
  if isfile(fno) == true
    rm(fno)
  end
  # --- Read data
  println(stderr," --- Read files")
  dat0, list0 = DelimitedFiles.readdlm(fn, header=true)
  dat1 = DelimitedFiles.readdlm(fn0)
  mean0 = Float64[]; std0 = Float64[]
  if all == false
    dat = dat0[1:end,1:6]
  else
    dat = copy(dat0)
  end
  mean0 = dat1[1:end,1]
  std0 = dat1[1:end,2]
  # --- Make list
  println(stderr," --- Make list")
  if all == false
    list = ["Displacement (EW) [m]","Displacement (NS) [m]","Displacement (UD) [m]","Shallow Gradient (EW) [sec/km]","Shallow Gradient (NS) [sec/km]","Gradient depth [km]"]
  else
    list = ["Displacement (EW) [m]","Displacement (NS) [m]","Displacement (UD) [m]","Shallow Gradient (EW) [sec/km]","Shallow Gradient (NS) [sec/km]","Gradient depth [km]","Long-term NTD 0_deg.","Long-term NTD 1_deg.","Long-term NTD 2_deg.","Long-term NTD 3_deg.","Long-term NTD 4_deg.","Scale factor 1","Scale factor 2"]
    for i in 1:NPB
      str0 = "3d B-spline NTD "
      str = string(str0,"$i [sec]")
      push!(list,str)
    end
  end
  # --- Plot
  println(stderr,"--- Drawing histograms")
  NN,MM = size(dat)
  dmin = minimum(dat,dims=1)
  dmax = maximum(dat,dims=1)
  dmean = mean(dat,dims=1)
  dmed = median(dat,dims=1)
  dstd = std(dat,dims=1)
  lmin = Float64[]; lmax = Float64[]
  if drawnls == false
    lmin = min.(dmin,dmean-5*dstd)
    lmax = max.(dmax,dmean+5*dstd)
  else
    for m in 1:MM
      if m<= 3
        push!(lmin,min(dmin[m],dmean[m]-5*dstd[m],mean0[m]-4*std0[m]))
        push!(lmax,max(dmax[m],dmean[m]+5*dstd[m],mean0[m]+4*std0[m]))
      elseif m >= 14
        push!(lmin,min(dmin[m],dmean[m]-5*dstd[m],mean0[m-10]-4*std0[m-10]))
        push!(lmax,max(dmax[m],dmean[m]+5*dstd[m],mean0[m-10]+4*std0[m-10]))
      else
        push!(lmin,min(dmin[m],dmean[m]-5*dstd[m]))
        push!(lmax,max(dmax[m],dmean[m]+5*dstd[m]))
      end
    end
  end
  for m in 1:MM
    println(stderr,"    Drawing histogram for $m / $MM")
    f(x) = 1/sqrt(2*pi*dstd[m]^2)*exp(-(x-dmean[m])^2/(2*dstd[m]^2))
    plt = histogram(dat[1:end,m],bins=range(dmin[m],dmax[m],length=nbins),normed = true,xlims=(lmin[m],lmax[m]),xlabel=list[m],ylabel="Density",lc=:white,fc=:skyblue,legend=:none,framestyle=:box,size=plot_size,left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin))
    plot!(plt,f,lw=2,lc=:blue)
    if drawnls == true
      mf = 0
      if m <= 3
        mf = m
      elseif m >= 14
        mf = m - 10
      end
      if m <= 3 || m >= 14
        g(x) = 1/sqrt(2*pi*std0[mf]^2)*exp(-(x-mean0[mf])^2/(2*std0[mf]^2))
        plot!(plt,g,lw=2,lc=:black)
      end
    end
    savefig(plt,"temp.pdf")
    PDFmerger.append_pdf!(fno,"temp.pdf",cleanup=true)
  end
  # --- Close process
  time2 = now()
  println(stderr," Start time:",time1)
  println(stderr," Finish time:",time2)
end
