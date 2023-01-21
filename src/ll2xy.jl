# If you change the default setting of GMT, you have to reset the setting especially for Projection
# gmt("gmtset PROJ_ELLIPSOID WGS-84")
#import GMT
#using Printf

export ll2xy, ll2xy_vec, xy2ll, xy2ll_vec, ll2xy_txt, xy2ll_txt
function ll2xy(lon,lat,lon0,lat0)
  deg1 = floor(lon0)
  deg2 = deg1 + 0.5
  x,y = GMT.gmt("mapproject -Jt$(lon0)/$(lat0)/1:1 -F -C -R$(deg1)/$(deg2)/-90/90",[lon lat])
  return x, y
end

function xy2ll(x,y,lon0,lat0)
  deg1 = floor(lon0)
  deg2 = deg1 + 0.5
  lon, lat = GMT.gmt("mapproject -Jt$(lon0)/$(lat0)/1:1 -I -F -C -R$(deg1)/$(deg2)/-90/90",[x y])
  return lon, lat
end

function ll2xy_vec(lon,lat,lon0,lat0)
  deg1 = floor(lon0)
  deg2 = deg1 + 0.5
  num = size(lon)[1]
  xy = GMT.gmt("mapproject -Jt$(lon0)/$(lat0)/1:1 -F -C -R$(deg1)/$(deg2)/-90/90",[lon lat])
  return xy[1:num,1],xy[1:num,2]
end

function xy2ll_vec(x,y,lon0,lat0)
  deg1 = floor(lon0)
  deg2 = deg1 + 0.5
  num = size(x)[1]
  ll = GMT.gmt("mapproject -Jt$(lon0)/$(lat0)/1:1 -I -F -C -R$(deg1)/$(deg2)/-90/90",[x y])
  return ll[1:num,1],ll[1:num,2]
end

function ll2xy_txt(fn0::String,fn::String,ks::Int64,lon0,lat0)
  deg1 = floor(lon0)
  deg2 = deg1 + 0.5
  f = open(fn0,"r")
  a = readlines(f)
  close(f)
  num1 = size(a)[1]
  num2 = size(split(a[1]))[1]
  a1 = Vector{String}(undef, num1)
  lon = zeros(num1)
  lat = zeros(num1)
  a2 = Vector{String}(undef, num1)
  if ks < 1
    error(" ll2xy_txt: ks must be larger than 1")
  elseif ks == 1
    if num2 == 2
      for i in 1:num1
        lon[i] = parse(Float64,split(a[i])[ks])
        lat[i] = parse(Float64,split(a[i])[ks+1])
      end
      x = GMT.gmt("mapproject -Jt$(lon0)/$(lat0)/1:1 -F -C -R$(deg1)/$(deg2)/-90/90",[lon lat])
      open(fn,"w") do out
        for i in 1:num1
          @printf(out,"%10.6f %10.6f\n",x[i,1],x[i,2])
        end
      end
    else
      for i in 1:num1
        lon[i] = parse(Float64,split(a[i])[ks])
        lat[i] = parse(Float64,split(a[i])[ks+1])
        a2[i] = join(split(a[i])[ks+2:num2]," ")
      end
      x = GMT.gmt("mapproject -Jt$(lon0)/$(lat0)/1:1 -F -C -R$(deg1)/$(deg2)/-90/90",[lon lat])
      open(fn,"w") do out
        for i in 1:num1
          @printf(out,"%10.6f %10.6f %s\n",x[i,1],x[i,2],a2[i])
        end
      end
    end
  elseif ks >= 2
    if ks+1 == num2
      for i in 1:num1
        a1[i] = join(split(a[i])[1:ks-1]," ")
        lon[i] = parse(Float64,split(a[i])[ks])
        lat[i] = parse(Float64,split(a[i])[ks+1])
     end
      x = GMT.gmt("mapproject -Jt$(lon0)/$(lat0)/1:1 -F -C -R$(deg1)/$(deg2)/-90/90",[lon lat])
      open(fn,"w") do out
        for i in 1:num1
          @printf(out,"%s %10.6f %10.6f\n",a1[i],x[i,1],x[i,2])
        end
      end
    else
      for i in 1:num1
        a1[i] = join(split(a[i])[1:ks-1]," ")
        lon[i] = parse(Float64,split(a[i])[ks])
        lat[i] = parse(Float64,split(a[i])[ks+1])
        a2[i] = join(split(a[i])[ks+2:num2]," ")
     end
      x = GMT.gmt("mapproject -Jt$(lon0)/$(lat0)/1:1 -F -C -R$(deg1)/$(deg2)/-90/90",[lon lat])
      open(fn,"w") do out
        for i in 1:num1
          @printf(out,"%s %10.6f %10.6f %s\n",a1[i],x[i,1],x[i,2],a2[i])
        end
      end
    end
  end
end 

function xy2ll_txt(fn0::String,fn::String,ks::Int64,lon0,lat0)
  deg1 = floor(lon0)
  deg2 = deg1 + 0.5
  f = open(fn0,"r")
  a = readlines(f)
  close(f)
  num1 = size(a)[1]
  num2 = size(split(a[1]))[1]
  a1 = Vector{String}(undef, num1)
  x = zeros(num1)
  y = zeros(num1)
  a2 = Vector{String}(undef, num1)
  if ks < 1
    error(" xy2ll_txt: ks must be larger than 1")
  elseif ks == 1
    if num2 == 2
      for i in 1:num1
        x[i] = parse(Float64,split(a[i])[ks])
        y[i] = parse(Float64,split(a[i])[ks+1])
      end
      ll = GMT.gmt("mapproject -Jt$(lon0)/$(lat0)/1:1 -I -F -C -R$(deg1)/$(deg2)/-90/90",[x y])
      open(fn,"w") do out
        for i in 1:num1
          @printf(out,"%3.9f %3.9f\n",ll[i,1],ll[i,2])
        end
      end
    else
      for i in 1:num1
        x[i] = parse(Float64,split(a[i])[ks])
        y[i] = parse(Float64,split(a[i])[ks+1])
        a2[i] = join(split(a[i])[ks+2:num2]," ")
      end
      ll = GMT.gmt("mapproject -Jt$(lon0)/$(lat0)/1:1 -I -F -C -R$(deg1)/$(deg2)/-90/90",[x y])
      open(fn,"w") do out
        for i in 1:num1
          @printf(out,"%3.9f %3.9f %s\n",ll[i,1],ll[i,2],a2[i])
        end
      end
    end
  elseif ks >= 2
    if ks+1 == num2
      for i in 1:num1
        a1[i] = join(split(a[i])[1:ks-1]," ")
        x[i] = parse(Float64,split(a[i])[ks])
        y[i] = parse(Float64,split(a[i])[ks+1])
     end
     ll = GMT.gmt("mapproject -Jt$(lon0)/$(lat0)/1:1 -I -F -C -R$(deg1)/$(deg2)/-90/90",[x y])
      open(fn,"w") do out
        for i in 1:num1
          @printf(out,"%s %3.9f %3.9f\n",a1[i],ll[i,1],ll[i,2])
        end
      end
    else
      for i in 1:num1
        a1[i] = join(split(a[i])[1:ks-1]," ")
        x[i] = parse(Float64,split(a[i])[ks])
        y[i] = parse(Float64,split(a[i])[ks+1])
        a2[i] = join(split(a[i])[ks+2:num2]," ")
      end
      ll = GMT.gmt("mapproject -Jt$(lon0)/$(lat0)/1:1 -I -F -C -R$(deg1)/$(deg2)/-90/90",[x y])
      open(fn,"w") do out
        for i in 1:num1
          @printf(out,"%s %3.9f %3.9f %s\n",a1[i],ll[i,1],ll[i,2],a2[i])
        end
      end
    end
  end
end 
