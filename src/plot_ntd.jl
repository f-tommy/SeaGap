#using Plots
#using DelimitedFiles

export plot_ntd
function plot_ntd(ntdrange=(-3,3),resrange=(-1,1); fno="ntd.pdf"::String,fn="residual.out"::String,autoscale=true,plot_size=(650,650),lmargin=2.5, rmargin=1.0, tmargin=1.0, bmargin=1.0, show=false,ms=4::Int64,bmargin0=-4)
  dat0 = DelimitedFiles.readdlm(fn)
  num = size(dat0)[1]
  dat = unixsort2(dat0,1,2)
  t0 = dat[1,1]
  numk = round(Int64,findmax(dat[1:num,2])[1])
  if autoscale == true
    p0 = plot(ylabel="NTD [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin0),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin),label="",xformatter=_->"")
    p1 = plot(xlabel="Time [hour]",ylabel="Residual [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, rmargin),legend=:none)
  else
    p0 = plot(ylim=ntdrange,ylabel="NTD [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin0),top_margin=Plots.Measures.Length(:mm, tmargin),right_margin=Plots.Measures.Length(:mm, rmargin),label="",xformatter=_->"")
    p1 = plot(ylim=resrange,xlabel="Time [hour]",ylabel="Residual [msec]",left_margin=Plots.Measures.Length(:mm, lmargin),bottom_margin=Plots.Measures.Length(:mm, bmargin),top_margin=Plots.Measures.Length(:mm, 0),right_margin=Plots.Measures.Length(:mm, rmargin),legend=:none)
  end
  for k in 1:numk
    t = Vector{Float64}(undef, 0)
    d = Vector{Float64}(undef, 0)
    dc = Vector{Float64}(undef, 0)
    dr = Vector{Float64}(undef, 0)
    i0 = 1
    for i in 1:num
      if round(Int64,dat[i,2]) == k
        push!(t,((dat[i,1] - t0) / (60*60)))
        push!(d,(dat[i,3] * 1000))
        push!(dc,(dat[i,4] * 1000))
        push!(dr,(dat[i,5] * 1000))
        i0 += 1
      end
    end
    scatter!(p0,t,d,ylabel="NTD [msec]",markershape=:cross,label="Transponder $k",markersize=ms)
    scatter!(p1,t,dr,xlabel="Time [hour]",ylabel="Residual [msec]",markershape=:cross,label="$k",markersize=ms)
  end
  t = (dat[1:num,1] .- t0) / (60*60)
  dc = dat[1:num,4] * 1000
  plot!(p0,t,dc,lc=:black,label="Smoothed NTD")
  plts = plot(p0,p1,layout=(2,1), size=plot_size,framestyle=:box)
  if show == false
    savefig(plts,fno)
  else
    gui(plts)
  end
end
