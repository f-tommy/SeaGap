#using Dates
#using Statistics
#using DelimitedFiles
#include("read_gnssa.jl")

function post_statistics(;fn0="obsdata.inp"::String,fn1="sample.out"::String,fn2="mcmc.out"::String,fno0="position.out"::String,fno1="statistics.out"::String,fno2="acceptance.out"::String)
  println(stderr," === Statistical post-processing for pos_array_mcmcpvg  ===")
  time1 = now()
  # --- Read data
  println(stderr," --- Read files")
  num, nk, tp, t1, x1, y1, z1, h1, p1, r1, t2, x2, y2, z2, h2, p2, r2, nf = read_obsdata(fn0)
  dat1, list = DelimitedFiles.readdlm(fn1,header=true)
  dat02 = DelimitedFiles.readdlm(fn2)
  dat2 = Int.(dat02[1:end,2:3])
  dat02 = dat2[dat2[:,1].==0,:][1:end,2]
  dat12 = dat2[dat2[:,1].==1,:][1:end,2]
  # --- Processing
  println(stderr," --- Processing")
  tave = mean(t1)
  ave1 = mean(dat1,dims=1)
  med1 = median(dat1,dims=1)
  std1 = std(dat1,dims=1)
  min1 = minimum(dat1,dims=1)
  max1 = maximum(dat1,dims=1)
  out = permutedims(vcat(list,ave1,std1,med1,min1,max1))
  acr0 = mean(dat02)
  acr1 = mean(dat12)
  # --- Output
  println(stderr," --- Output")
  open(fno0,"w") do out0
    println(out0,"$tave $(ave1[1]) $(ave1[2]) $(ave1[3]) $(std1[1]) $(std1[2]) $(std1[3])")
  end
  open(fno1,"w") do out1
    println(out1,"#Parameter mean std median min max")
    Base.print_array(out1,out)
    println(out1,"")
  end
  open(fno2,"w") do out2
    println(out2,"Acceptance_ratio_MCMC-1: $acr0")
    println(out2,"Acceptance_ratio_MCMC-2: $acr1")
  end
  # --- Close process
  time2 = now()
  println(stderr," Start time:",time1)
  println(stderr," Finish time:",time2)
  println(stderr," Check $fno1 $fno2")
end
