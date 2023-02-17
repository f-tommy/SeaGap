#=
  Functions for reading GNSS-A data
  written by Fumiaki Tomita 2022/06/29
=#
#using DelimitedFiles

export read_matrix,read_info, read_prof, read_ant, read_gps, read_initial, read_jttq, read_pxppos, read_obsdata
# === Read matrix
"""
    read_matrix(filename)

 Read text data file; then, return its size (`N`*`M`) and matrix (`a`)

# Example
    N, M, a = read_matrix("site_info.txt")

"""
function read_matrix(filename::String)
  a = DelimitedFiles.readdlm(filename)
  n, m = size(a)
  println(stderr," --- Read $filename: $n")
  return n,m,a
end 

# === Read site info
"""
    read_info(site,filename)

Read site location data `filename` and return the location for a certain site given by `site`
`filename` is a text file with (1:Site_name(str), 2:Lon(deg), 3:Lat(deg), 4:Water_depth(m), 5:Total number of transponders(Int)).

# Example
    lon,lat,dep,numk = read_info("G20","site_info.txt")
"""
function read_info(site::String,filename="site_info.txt")
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
"""
    read_prof(filename,TR_DEPTH)

Read a sound speed profile `filename`.

* `filename`: Input file name for a sound speed profile
* `TR_DEPTH`: Transducer depth from the sea-surface

Output:
* `z`: Arrangement for depth
* `v`: Arrangement for velocity
* `nz_st`: Layer number which inculdes transducer
* `numz`: Number of layers

# Example
    z, v, nz_st, numz = read_prof("ss_prof.zv",3.0)
"""
function read_prof(filename::String,TR_DEPTH)
  # --- Read sound speed data
  a = DelimitedFiles.readdlm(filename)
  numz = size(a)[1]
  z = a[1:numz,1]
  v = a[1:numz,2]
  # --- Find nz_st
  nz = 1
  nz_st = 1
  while z[nz] < TR_DEPTH
    nz_st = nz
    nz += 1
  end
  println(" --- Read $filename: $nz_st, $numz")
  return z,v,nz_st,numz
end

# === Read TR-ANT offset
"""
    read_ant(filename)

Read a file `filename` for an offset between GNSS antenna and transducer, and store it as an arrangement `e`.

# Example
    e = read_ant("tr-ant.inp")
"""
function read_ant(filename::String)
  e0 = DelimitedFiles.readdlm(filename)
  e = transpose(e0)
  println(stderr," --- Read $filename: $(e[1]), $(e[2]), $(e[3])")
  return e
end

# === Read gps positions
"""
    read_gps(filename)

Read a file `filename` for time-series of GNSS antenna positions and attitudes, and store them as arrangements.

Output:
* `numj`: Total number of data
* `tg0`: Time [sec]
* `xg0`: EW GNSS positon [m]
* `yg0`: NS GNSS positon [m]
* `zg0`: UD GNSS positon [m]
* `hd`: Heading [deg]
* `pd`: Pitch [deg]
* `rd`: Roll [deg]

# Example
    numj,tg0,xg0,yg0,zg0,hd,pd,rd = read_gps("gps.jxyhhpr")
"""
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

"""
    read_initial(filename)

Read a file `filename` for initial values used in `static_array_mcmcgrad` or `static_array_mcmcgradc`, and store them as arrangements.

Output:
* `np`: Number of data
* `a0`: Initial value
* `a1`: Lower limit 
* `a2`: Upper limit 
* `da`: Step width
* `list`: Parameter name

# Example
    np, a0, a1, a2, da, list = read_initial("initial.inp")
"""
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
"""
    read_jttq(filename,k,tp,tt0,tt,qq)

Read a file `filename` of travel-time data for `k`th transponder, and store them as the arrangements `tp`, `tt0`, `tt`, and `qq`.

Output:
* `nump`: Total number of acoustic shots
* `tp`: Time at each shot [sec]
* `tt0`: Travel-time with mechanical delay
* `tt`: Travel-time without delay
* `qq`: Quality value (this is not used; thus, you can put zero for all, if you do not have any quality identifers)

# Example
    maxp = 5000
    nump = zeos(3); tp = zeros(3,maxp); tt0 = zeros(3,maxp)
    tt = zeros(3,maxp); qq = zeros(3,maxp)
    nump[1] = read_jttq("pxp-1.jttq",1,tp,tt0,tt,qq)
    nump[2] = read_jttq("pxp-2.jttq",2,tp,tt0,tt,qq)
    nump[3] = read_jttq("pxp-3.jttq",3,tp,tt0,tt,qq)
"""
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
"""
    read_pxppos(filename)

Read a file `filename` for seafloor transponder positions, and store them as arrangements.

Output:
* `numk`: Total number of the seafloor transponders
* `px`: `numk`-vector of EW postion [m]
* `py`: `numk`-vector of NS postion [m]
* `pz`: `numk`-vector of UD postion [m]

# Example
    numk, px, py ,pz = read_pxppos("pxp-ini.xyh")
"""
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
"""
    read_obsdata(filename)

Read a file `filename` for seafloor transponder positions, and store them as arrangements.

Output:
* `num`: Total number of data
* `nk`: Transponder number
* `tp`: Travel-Time [sec]
* `t1`: Signal transmitting time [sec]
* `x1`: EW GNSS antenna position when transmitting [m]
* `y1`: NS GNSS antenna position when transmitting [m]
* `z1`: UD GNSS antenna position when transmitting [m]
* `h1`: Heading when transmitting [deg]
* `p1`: Pitch when transmitting [deg]
* `r1`: Roll when transmitting [deg]
* `t2`: Signal recieving time [sec]
* `x2`: EW GNSS antenna position when recieving [m]
* `y2`: NS GNSS antenna position when recieving [m]
* `z2`: UD GNSS antenna position when recieving [m]
* `h2`: Heading when recieving [deg]
* `p2`: Pitch when recieving [deg]
* `r2`: Roll when recieving [deg]
* `nf`: Shot group number for `kinematic_array()`

# Example
    num,nk,tp,t1,x1,y1,z1,h1,p1,r1,t2,x2,y2,z2,h2,p2,r2,nf = read_obsdata("obsdata.inp")
"""
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
  nf = Int.(round.(a[1:num,17])) # Flag for kinematic_array
  println(stderr," --- Read $filename: $num")
  if num < 3
    error(" read_obsdata: number of lines must be more than 3")
  end
  return num,nk,tp,t1,x1,y1,z1,h1,p1,r1,t2,x2,y2,z2,h2,p2,r2,nf
end

"""
    read_ntd(filename)

Read a file `filename` for NTD estimation results obtained by `static_array()` for example, and store them as arrangements.

Output:
* `ts`: Shot Time [sec]
* `nk`: Transponder number
* `to`: Projected travel-time residuals in the nadir direction
* `tc`: NTD modeled by B-spline bases
* `tr`: Travel-time residulas subtracting `tc` from `to` 

# Example
    ts, nk, to, tc, tr = read_ntd("ntd.out")
"""
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
