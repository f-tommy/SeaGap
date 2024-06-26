#using Plots
#using DelimitedFiles

export plot_bspline_gradv
function plot_bspline_gradv(mode=1::Int64,lr=(-0.2,0.2),sr=(-0.2,0.2),sr1=(0.0,0.1),sr2=(0.0,0.1),dr1=(0.0,2.0),dr2=(0.0,2.0),gr1=(0.0,0.1),gr2=(0.0,0.1);fn0="synthetic_NTD.txt"::String,fn="bspline_grad.out"::String,fno="bspline_all.pdf"::String,plot_size=(1000,1000),lmargin=4.0,tmargin=1.0,bmargin=0.5,rmargin=1.0,bmargin0=-3.0,show=false::Bool,ms=5::Int64,gfs=10::Int64,leg="Estimated"::String,autoscale=true::Bool)
  if mode == 1
    n, m, dat1 = read_matrix(fn)
    t1 = (dat1[:,1].-dat1[1,1])/3600
    lntd1 = dat1[:,2]*1000
    sntd1 = dat1[:,3]*1000
    gsew1 = dat1[:,4]*1000
    gsns1 = dat1[:,5]*1000
    gdew1 = dat1[:,6]*1
    gdns1 = dat1[:,7]*1
    ggew1 = gdew1.*gsew1./2.0
    ggns1 = gdns1.*gsns1./2.0
  elseif mode == 2
    n, m, dat1 = read_matrix(fn0)
    t1 = (dat1[:,1].-dat1[1,1])/3600
    lntd1 = dat1[:,2]*1000
    sntd1 = dat1[:,3]*1000
    gsew1 = dat1[:,4]*1000
    gsns1 = dat1[:,5]*1000
    gdew1 = dat1[:,6]*1
    gdns1 = dat1[:,7]*1
    ggew1 = gdew1.*gsew1./2.0
    ggns1 = gdns1.*gsns1./2.0
    leg = "Synthetic"
  else
    n, m, dat1 = read_matrix(fn)
    t1 = (dat1[:,1].-dat1[1,1])/3600
    lntd1 = dat1[:,2]*1000
    sntd1 = dat1[:,3]*1000
    gsew1 = dat1[:,4]*1000
    gsns1 = dat1[:,5]*1000
    gdew1 = dat1[:,6]*1
    gdns1 = dat1[:,7]*1
    ggew1 = gdew1.*gsew1./2.0
    ggns1 = gdns1.*gsns1./2.0
    n, m, dat2 = read_matrix(fn0)
    t2 = (dat2[:,1].-dat2[1,1])/3600
    lntd2 = dat2[:,2]*1000
    sntd2 = dat2[:,3]*1000
    gsew2 = dat2[:,4]*1000
    gsns2 = dat2[:,5]*1000
    gdew2 = dat2[:,6]*1
    gdns2 = dat2[:,7]*1
    ggew2 = gdew2.*gsew2./2.0
    ggns2 = gdns2.*gsns2./2.0
  end
  if autoscale == false
    p1 = plot(t1,lntd1,ylims=lr,xlabel="Time [h]",label=leg,ylabel="L-NTD [ms]",framestyle=:box)
    p2 = plot(t1,sntd1,ylims=sr,xlabel="Time [h]",label=leg,ylabel="S-NTD [ms]",framestyle=:box)
    p3 = plot(t1,gsew1,ylims=sr1,xlabel="Time [h]",label=leg,ylabel="S-Grad (EW) [ms/km]",framestyle=:box)
    p4 = plot(t1,gsns1,ylims=sr2,xlabel="Time [h]",label=leg,ylabel="S-Grad (NS) [ms/km]",framestyle=:box)
    p5 = plot(t1,gdew1,ylims=dr1,xlabel="Time [h]",label=leg,ylabel="Grad Depth (EW) [km]",framestyle=:box)
    p6 = plot(t1,gdns1,ylims=dr2,xlabel="Time [h]",label=leg,ylabel="Grad Depth (NS) [km]",framestyle=:box)
    p7 = plot(t1,ggew1,ylims=gr1,xlabel="Time [h]",label=leg,ylabel="D-Grad (EW) [ms]",framestyle=:box)
    p8 = plot(t1,ggns1,ylims=gr2,xlabel="Time [h]",label=leg,ylabel="D-Grad (NS) [ms]",framestyle=:box)
  else
    p1 = plot(t1,lntd1,xlabel="Time [h]",label=leg,ylabel="L-NTD [ms]",framestyle=:box)
    p2 = plot(t1,sntd1,xlabel="Time [h]",label=leg,ylabel="S-NTD [ms]",framestyle=:box)
    p3 = plot(t1,gsew1,xlabel="Time [h]",label=leg,ylabel="S-Grad (EW) [ms/km]",framestyle=:box)
    p4 = plot(t1,gsns1,xlabel="Time [h]",label=leg,ylabel="S-Grad (NS) [ms/km]",framestyle=:box)
    p5 = plot(t1,gdew1,xlabel="Time [h]",label=leg,ylabel="Grad Depth (EW) [km]",framestyle=:box)
    p6 = plot(t1,gdns1,xlabel="Time [h]",label=leg,ylabel="Grad Depth (NS) [km]",framestyle=:box)
    p7 = plot(t1,ggew1,xlabel="Time [h]",label=leg,ylabel="D-Grad (EW) [ms]",framestyle=:box)
    p8 = plot(t1,ggns1,xlabel="Time [h]",label=leg,ylabel="D-Grad (NS) [ms]",framestyle=:box)
  end
  if mode == 3
    plot!(p1,t2,lntd2,label="Synthetic")
    plot!(p2,t2,sntd2,label="Synthetic")
    plot!(p3,t2,gsew2,label="Synthetic")
    plot!(p4,t2,gsns2,label="Synthetic")
    plot!(p5,t2,gdew2,label="Synthetic")
    plot!(p6,t2,gdns2,label="Synthetic")
    plot!(p7,t2,ggew2,label="Synthetic")
    plot!(p8,t2,ggns2,label="Synthetic")
  end
  plt = plot(p1,p2,p3,p4,p5,p6,p7,p8,layout=(4,2),size=plot_size)
#  plt = plot(p1,p3,p5,p7,p2,p4,p6,p8,layout=(2,4),size=plot_size)
  if show == false
    savefig(plt,fno)
  else
    gui(plt)
  end
end
