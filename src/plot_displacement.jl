#using Plots
#using DelimitedFiles
#using Statistics
#using Printf
#include("/usr/local/share/julia_bin/unixsort.jl")
#include("LineFitting.jl")

export plot_displacement

"""
    plot_displacement(ts,te,EW_range,NS_range,UD_range; autoscale,sigma,cal,weight,predict,fno1,fno2,fno3,fn,plot_size,lmargin,rmargin,tmargin,bmargin,bmargin0,show,ms,lfs,tfs,int,lw,msw,pscale,fa,lgs)

Make a figure of array displacement time-series.
The input file `fn` is converted time-series data through `convert_displacement()`. 
Optionally, you can calculate linear fitting lines.

* `ts`: Start time of the figure (Time at the left edge of the figure) [yr]
* `te`: End time of the figure (Time at the right edge of the figure) [yr]
* `EW_range`, `NS_range`, `UD_range`: If `autoscale=false`, the range of Y-axis is set to be those ranges
* autoscale: if `autoscale=true` (default), the range of Y-axis is automatically determined depending on the usable array displacements
* `sigma`: Errorbar scaling (Errorbar scaling (Since standard deviation of an array displacement tends to be small, the plotting errorbar is multiplied by `sigma`: `sigma=10` by default)
* `cal`: if `cal=true` (default), a regression line is calculated for each component, is shown in the figure, and the regression results are written in `fno2` as a text file
* `weight`: if `weight=true` (default), the regression line is calculated considering the weight; the weight is provided as an inverse square of the observation error (`weight=true` by default)
* `predict`: if `predict=true` (default), the predicted values from the regression line are written in "`fno3`-EW.txt", "`fno3`-NS.txt", and "`fno3`-UD.txt"
* `fno1`: Output figure name (`fno1="time-series.pdf"` by default)
* `fno2`: Name of an output text file showing array displacement rates (`fno2="velocity.out"` by default)
* `fno3`: Indicator name of output text files showing the predicted values from the regression lines (`fno3="predict"` by default)
* `fn`: Input file name (`fn="converted_position.out"` by default)
* `show`:  If `show=false`, the figure is saved as `fno` (`fno` is name of the output figure). If `show=true` in REPL, a figure is temporally shown.
* `plot_size`: Figure size (`plot_size=(550,650)` by default)
* `lmargin`: Plot margin for the left edge (`lmargin=1.5` by default)
* `rmargin`: Plot margin for the right edge (`rmargin=1.5` by default)
* `tmargin`: Plot margin for the top edge (`tmargin=1.0` by default)
* `bmargin`: Plot margin for the bottom edge (`bmargin=-1.0` by default)
* `bmargin0`: Plot margin for the bottom edges of upper two panels (`bmargin0=-4.5` by default)
* `show`: if `show=true`, a figure is temporally shown; if false, the figure is saved as `fno` (`show=false` by default)
* `ms`: Plotted marker size (`ms=3` by default)
* `lfs`: Fontsize for label (`lfs=9` by default)
* `tfs`: Fontsize for tick (`tfs=7` by default)
* `lw`: Width of fitting lines (`lw=2` by default)
* `msw`: Width of marker outline (Errorbar thickness, `msw=1` by default)
* `pscale`: if you set smaller value, the autoscale Y-axis range is wider (`pscale=0.5` by default)
* `fa`: Transparency of fill color of confidence intervali (`fa=0.5` by default)
* `lgs`: Label Fontsize (`lgs=6` by default)
* `alp`: Percentage for confidence interval (`alp=0.95` by default)

# Example
    plot_displacement(2012,2019,pscale=0.5)
"""
function plot_displacement(ts,te,EW_range=(-1,1),NS_range=(-1,1),UD_range=(-1,1); autoscale=true,sigma=10,cal=true::Bool,weight=true::Bool,predict=true::Bool,fno1="time-series.pdf"::String,fno2="velocity.out",fno3="predict",fn="converted_position.out"::String,plot_size=(550,650),lmargin=1.5,rmargin=1.5, tmargin=1.0, bmargin=-1,bmargin0=-4.5,show=false::Bool,ms=3::Int64,lfs=9::Int64,tfs=7::Int64,lw=2,msw=1,pscale=0.5,fa=0.25,lgs=6::Int64,alp=0.95)
  # --- Read
  println(stderr," --- Read $fn")
  dat0, header = DelimitedFiles.readdlm(fn,header=true)
  dat = unixsort(dat0[:,1:10],1)
  num = size(dat)[1]
  t = dat[:,1]; dx0 = dat[:,2]; dy0 = dat[:,3]; dz0 = dat[:,4]; sx0 = dat[:,5]; sy0 = dat[:,6]; sz0 = dat[:,7]
  idx = Int.(floor.(dat[:,8])); idy = Int.(floor.(dat[:,9])); idz = Int.(floor.(dat[:,10]))
  tx = t[idx.==1]; ty = t[idy.==1]; tz = t[idz.==1]
  dx = dx0[idx.==1]; dy = dy0[idy.==1]; dz = dz0[idz.==1]
  sx = sx0[idx.==1]; sy = sy0[idy.==1]; sz = sz0[idz.==1]
  nx = length(tx); ny = length(ty); nz = length(tz)
  println(stderr,"  Number of data, X:$nx Y:$ny Z:$nz")
  # --- Calculate velocity
  if cal == true
    println(stderr," --- Calculate velocity")
    if nx >= 3 && ny >= 3 && nz >= 3
      newTx = collect(range(tx[1],tx[end],50))
      newTy = collect(range(ty[1],ty[end],50))
      newTz = collect(range(tz[1],tz[end],50))
      if weight == true
        lsx = linefit(tx,dx,sx.^(-2),alpha=alp,newX=newTx)
        lsy = linefit(ty,dy,sy.^(-2),alpha=alp,newX=newTy)
        lsz = linefit(tz,dz,sz.^(-2),alpha=alp,newX=newTz)
      else
        lsx = linefit(tx,dx,alpha=alp,newX=newTx)
        lsy = linefit(ty,dy,alpha=alp,newX=newTy)
        lsz = linefit(tz,dz,alpha=alp,newX=newTz)
      end
      vx = lsx.coef[2]; vy = lsy.coef[2]; vz = lsz.coef[2]
      svx = lsx.coefstd[2]; svy = lsy.coefstd[2]; svz = lsz.coefstd[2]
      rmsx = sqrt(lsx.sigma2); rmsy = sqrt(lsy.sigma2); rmsz = sqrt(lsz.sigma2)
      cxy = cor(lsx.misfit,lsy.misfit)
      velx = string(round(vx,digits=1),"\u00B1",round(svx,digits=1)," cm/yr")
      vely = string(round(vy,digits=1),"\u00B1",round(svy,digits=1)," cm/yr")
      velz = string(round(vz,digits=1),"\u00B1",round(svz,digits=1)," cm/yr")
      open(fno2,"w") do out2
        println(out2,"# (Vx,Vy,Vz,Sx,Sy,Sz,RMSx,RMSy,RMSz,CorXY) in cm Linear model")
        println(out2,"$vx $vy $vz $svx $svy $svz $rmsx $rmsy $rmsz $cxy $nx $ny $nz")
      end
      if predict == true
        fno3x = string(fno3,"-EW.txt")
        fno3y = string(fno3,"-NS.txt")
        fno3z = string(fno3,"-UD.txt")
        open(fno3x,"w") do out3x
          println(out3x,"# Year Disp. Pred. Residual Sigma")
          Base.print_array(out3x,hcat(tx,dx,lsx.pred,lsx.misfit,sx))
          println(out3x,"")
        end
        open(fno3y,"w") do out3y
          println(out3y,"# Year Disp. Pred. Residual Sigma")
          Base.print_array(out3y,hcat(ty,dy,lsy.pred,lsy.misfit,sy))
          println(out3y,"")
        end
        open(fno3z,"w") do out3z
          println(out3z,"# Year Disp. Pred. Residual Sigma")
          Base.print_array(out3z,hcat(tz,dz,lsz.pred,lsz.misfit,sz))
          println(out3z,"")
        end
      end
      println(stderr," Output velocity: $fno2")
    else
      error("Veclocity cannot be calculated: number of data $nx $ny $nz should be >= 3")
    end
  end
  # --- Plot
  println(stderr," --- Plot")
  minx = minimum(dx); miny = minimum(dy); minz = minimum(dz)
  maxx = maximum(dx); maxy = maximum(dy); maxz = maximum(dz)
  deltax = maxx - minx; deltay = maxy - miny; deltaz = maxz - minz
  lowx = minx - deltax/pscale; uppx = maxx + deltax/pscale
  lowy = miny - deltay/pscale; uppy = maxy + deltay/pscale
  lowz = minz - deltaz/pscale; uppz = maxz + deltaz/pscale
  if autoscale == true
    pltx = plot(ylims=(lowx,uppx)); plty = plot(ylims=(lowy,uppy)); pltz = plot(ylims=(lowz,uppz))
  else
    pltx = plot(ylims=EW_range);  plty = plot(ylims=NS_range); pltz = plot(ylims=UD_range)
  end
  if cal == true
    plot!(pltx,newTx,lsx.pred_new,ribbon=(lsx.pred_new - lsx.pred_lower,lsx.pred_upper - lsx.pred_new),label=velx,lc=:skyblue,lw=lw,fillalpha=fa)
    plot!(plty,newTy,lsy.pred_new,ribbon=(lsy.pred_new - lsy.pred_lower,lsy.pred_upper - lsy.pred_new),label=vely,lc=:skyblue,lw=lw,fillalpha=fa)
    plot!(pltz,newTz,lsz.pred_new,ribbon=(lsz.pred_new - lsz.pred_lower,lsz.pred_upper - lsz.pred_new),label=velz,lc=:skyblue,lw=lw,fillalpha=fa)
  end
  scatter!(pltx,t,dx0,ylabel="Disp. (EW) [cm]",label="",xformatter=_->"",markershape=:square,xlims=(ts,te),framestyle=:box,ms=ms,mc=RGB(0.7,0.7,0.7),wsc=:black,top_margin=Plots.Measures.Length(:mm,tmargin),bottom_margin=Plots.Measures.Length(:mm,bmargin0),msw=0.5)
  scatter!(plty,t,dy0,ylabel="Disp. (NS) [cm]",xformatter=_->"",label="",markershape=:square,xlims=(ts,te),framestyle=:box,ms=ms,mc=RGB(0.7,0.7,0.7),wsc=:black,top_margin=Plots.Measures.Length(:mm,0),bottom_margin=Plots.Measures.Length(:mm,bmargin0),msw=0.5)
  scatter!(pltz,t,dz0,xlabel="Time [yr]",ylabel="Disp. (UD) [cm]",label="",markershape=:square,xlims=(ts,te),framestyle=:box,ms=ms,mc=RGB(0.7,0.7,0.7),wsc=:black,top_margin=Plots.Measures.Length(:mm,0),bottom_margin=Plots.Measures.Length(:mm,bmargin),msw=0.5)
  scatter!(pltx,tx,dx,yerror=sx*sigma,label="",markershape=:square,ms=ms,mc=:black,wsc=:black,msw=msw)
  scatter!(plty,ty,dy,yerror=sy*sigma,label="",markershape=:square,ms=ms,mc=:black,wsc=:black,msw=msw)
  scatter!(pltz,tz,dz,yerror=sz*sigma,label="",markershape=:square,ms=ms,mc=:black,wsc=:black,msw=msw)
  plt = plot(pltx,plty,pltz,layout=(3,1),size=plot_size,link=:both,labelfontsize=lfs,tickfontsize=tfs,left_margin=Plots.Measures.Length(:mm,lmargin),right_margin=Plots.Measures.Length(:mm, rmargin),legendfont=font(lgs,"Times"))
  println(stderr," --- Save")
  if show == false
    savefig(plt,fno1)
  else
    gui(plt)
  end
end
