#=
  Functions for reading GNSS-A data
  written by Fumiaki Tomita 2022/06/29
=#
#using DelimitedFiles

export read_info, read_prof, read_ant, read_gps, read_initial, read_jttq, read_pxppos, read_obsdata
# === Read site info
function read_info(filename::String,site::String)
  # 1:Site_name(str), 2:Lon(deg), 3:Lat(deg), 4:Water_depth(m), 5:total_number_of_transponders(Int)
  f = open(filename,"r")
  a = readlines(f)
  close(f)
  num = size(a)[1]
  for i in 1:num
    b = split(a[i])
    if b[1] == site
      lon = parse(Float64,b[2])
      lat = parse(Float64,b[3])
      dep = parse(Float64,b[4])
      numk = parse(Int64,b[5])
      return lon,lat,dep,numk
    end
  end
end
    
# === Read sound speed file
function read_prof(filename::String,XDUCER_DEPTH)
  # --- Read sound speed data
  a = DelimitedFiles.readdlm(filename)
  numz = size(a)[1]
  z = a[1:numz,1]
  v = a[1:numz,2]
  # --- Find nz_st
  nz = 1
  nz_st = 1
  while z[nz] < XDUCER_DEPTH
    nz_st = nz
    nz += 1
  end
  println(" --- Read $filename: $nz_st, $numz")
  return z,v,nz_st,numz
end

# === Read TR-ANT offset
function read_ant(filename::String)
  e0 = DelimitedFiles.readdlm(filename)
  e = transpose(e0)
  println(stderr," --- Read $filename: $(e[1]), $(e[2]), $(e[3])")
  return e
end

# === Read gps positions
function read_gps(filename::String)
  a = DelimitedFiles.readdlm(filename)
  numj = size(a)[1]
  tg0 = a[1:numj,1]
  xg0 = a[1:numj,2]
  yg0 = a[1:numj,3]
  zg0 = a[1:numj,4]
  hd = a[1:numj,5]
  pd = a[1:numj,6]
  rd = a[1:numj,7]
  println(stderr," --- Read $filename: $numj")
  return numj,tg0,xg0,yg0,zg0,hd,pd,rd
end

function read_initial(filename::String)
  a = DelimitedFiles.readdlm(filename)
  np = size(a)[1]
  a0 = Float64.(a[1:np,1])
  a1 = Float64.(a[1:np,2])
  a2 = Float64.(a[1:np,3])
  da = Float64.(a[1:np,4])
  list = a[1:np,5]
  println(stderr," --- Read $filename: $np")
  if np < 14
    error("  read_initial: number of lines must be more than 14")
  end
  return np, a0, a1, a2, da, list
end

# === Read jttq positions
function read_jttq(filename::String,k::Int64,tp,tt0,tt,qq)
  a = DelimitedFiles.readdlm(filename)
  nump0 = size(a)[1]
  tp[k,1:nump0] = a[1:nump0,1]
  tt0[k,1:nump0] = a[1:nump0,2]
  tt[k,1:nump0] = a[1:nump0,3]
  qq[k,1:nump0] = a[1:nump0,4]
  println(stderr," --- Read $filename: $nump0")
  return nump0
end

# === Read PXP positions
function read_pxppos(filename::String)
  a = DelimitedFiles.readdlm(filename)
  numk = size(a)[1]
  px = a[1:numk,1]
  py = a[1:numk,2]
  pz = a[1:numk,3]
  println(stderr," --- Read $filename: $numk")
  for k in 1:numk
    println(stderr,"          PxP $k: $(px[k]), $(py[k]), $(pz[k])")
  end
  return numk,px,py,pz
end 

# === Read observational file as described by JCG-style
function read_obsdata(filename::String)
  # --- Read observation data
  a = DelimitedFiles.readdlm(filename)
  num = size(a)[1]
  nk = Int.(round.(a[1:num,1]))  # PxP num
  tp = a[1:num,2]                # Obs TT
  t1 = a[1:num,3]                # Shot time
  x1 = a[1:num,4]                # Shot pos-X
  y1 = a[1:num,5]                # Shot pos-Y
  z1 = a[1:num,6]                # Shot pos-Z
  h1 = a[1:num,7]                # Shot heading
  p1 = a[1:num,8]                # Shot pitch
  r1 = a[1:num,9]                # Shot roll
  t2 = a[1:num,10]               # Recieve time
  x2 = a[1:num,11]               # Recieve pos-X
  y2 = a[1:num,12]               # Recieve pos-Y
  z2 = a[1:num,13]               # Recieve pos-Z
  h2 = a[1:num,14]               # Recieve heading
  p2 = a[1:num,15]               # Recieve pitch
  r2 = a[1:num,16]               # Receive roll
  nf = Int.(round.(a[1:num,17])) # Flag for pos_array_each
  println(stderr," --- Read $filename: $num")
  if num < 3
    error(" read_obsdata: number of lines must be more than 3")
  end
  return num,nk,tp,t1,x1,y1,z1,h1,p1,r1,t2,x2,y2,z2,h2,p2,r2,nf
end

function read_ntd(filename::String)
  # --- Read observation data
  a = DelimitedFiles.readdlm(filename)
  ts = a[1:end,1]                # Shot Time
  nk = Int.(round.(a[1:end,2]))  # PxP num
  to = a[1:end,3]                # NTD
  tc = a[1:end,4]                # Spline-Modeled NTD
  tr = a[1:end,5]                # Residual
  println(stderr," --- Read $filename")
  return ts, nk, to, tc, tr
end
