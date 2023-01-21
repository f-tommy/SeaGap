#=
  Functions for travel-time calculation 
    localradius : function for calculating Local Earth Radius (Chadwell & Sweeney, 2010)
    ttvert0 : function for calculating vertical travel-time ()
    ttsphere : function for calculating simple travel-time (See Tomita & Kido, in prep)
    ttsphere : function for calculating simple travel-time (See Tomita & Kido, in prep)
    ttcorrection : function to obtain correction terms for approximate travel-time (See Tomita & Kido, in prep)
    xyz2tt : function for exact travel-time calculation
    xyz2tt_rapid : function for approximate travel-time calculation
  written by Fumiaki Tomita 2022/06/29
=#
#using DelimitedFiles
#using LinearAlgebra
#using Statistics
export localradius, ttvert0, ttsphere, xyz2tt, ttcorrection, xyz2tt_rapid, xyz2ttg_rapid

# === Local Earth radius
function localradius(lat)
  a = 6378137.0 # WGS84 long radius
  b = 6356752.31424 # WGS84 short radius
  e2 = (a^2.0-b^2.0) / a^2.0
  rlat2 = sin(lat/180.0*pi)^2.0
  Rm = a*(1.0-e2)/sqrt((1.0-e2*rlat2)^3.0)
  Rn = a/sqrt(1.0-e2*rlat2)
  Rg = sqrt(Rm*Rn) # Gaussian mean Earth radius
  Rl = 1.0/ ((1.0/Rm+1.0/Rn)/2.0)  # Local mean Earth radius
  return Rg, Rl
end

# === Calculte tt_vert0
function ttvert0(pxp_height,xducer_height,XDUCER_DEPTH,z,v,nz_st,numz::Int64)
  v_range = xducer_height - pxp_height
  pxp_depth = XDUCER_DEPTH + v_range
  tt_vert0 = 0.0
  id = 0
  v_deep = 0.0
  nz = nz_st
  while id == 0
    if nz == nz_st
      delta_z = (z[nz+1] - XDUCER_DEPTH)
    else
      delta_z = (z[nz+1] - z[nz])
    end
    if z[nz+1] >= pxp_depth
      id = 1
      delta_z = (pxp_depth - z[nz])
      grad = (v[nz+1] - v[nz]) / (z[nz+1] - z[nz])
      v_deep = v[nz] + grad*delta_z
    end
    tt_vert0 += delta_z / v[nz]
    nz += 1
  end
  return tt_vert0, v_deep, v_range
end

# === Calculate tt_sphere
function ttsphere(pxpx,pxpy,pxpz,rsx,rsy,rsz,Ra,tt_vert0,v_deep,v_range)
  delta_z = rsz - pxpz
  theta = sqrt((rsx-pxpx)^2.0+(rsy-pxpy)^2.0) / Ra
  b = Ra + rsz
  c = Ra + pxpz
  a = sqrt(b^2.0+c^2.0-2.0*b*c*cos(theta))
  delta_h = delta_z - v_range
  tt_vert = tt_vert0 + delta_h / v_deep
  tt_sphere = tt_vert*a/delta_z
  return tt_sphere, theta
end

# === Calculate traveltime
function xyz2tt(pxpx,pxpy,pxpz,rsx,rsy,rsz,z,v,nz_st::Int64,numz::Int64,Ra,XDUCER_DEPTH)
  # --- Set basic parameters
  eps_dist = 1.e-5
  ITMAX = 50
  # --- Set parameter
  xducer_height = rsz
  pxp_height = pxpz
  pxp_depth = XDUCER_DEPTH + xducer_height - pxp_height
  b = Ra + xducer_height
  c = Ra + pxp_height
  x = sqrt((pxpx-rsx)^2+(pxpy-rsy)^2)
  theta = x/Ra
  s_range = sqrt(b^2+c^2-2.0*b*c*cos(theta))
  vert = (rsz-pxpz)/s_range
  sinxi0 = c/s_range*sin(theta)
  p = b/v[1]*sinxi0
  # --- Loop calculation
  n = 0
  tt = 0.e0 # Essential to define local variable
  diff_dist = 1.e6
  while diff_dist >= eps_dist
    if n > ITMAX
      break
    end
    n += 1
    tt = 0.e0
    delta = 0.e0
    ddelta = 0.e0
    nz = nz_st
    id = 0
    while id == 0
      if nz == nz_st
        r0 = Ra + xducer_height
        r1 = Ra + xducer_height - z[nz+1] + XDUCER_DEPTH
      else
        r0 = Ra + xducer_height + XDUCER_DEPTH - z[nz]
        r1 = Ra + xducer_height + XDUCER_DEPTH - z[nz+1]
      end
      if z[nz+1] >= pxp_depth
        id = 1
        r1 = Ra + pxp_height
      end
      delta += acos(v[nz]*p/r0) - acos(v[nz]*p/r1)
      ddelta += -1.0/sqrt((r0/v[nz])^2.0-p^2.0) + 1.0/sqrt((r1/v[nz])^2.0-p^2.0)
      tt += sqrt((r0/v[nz])^2.0-p^2.0) - sqrt((r1/v[nz])^2.0-p^2.0)
      nz += 1
    end
    diff_theta = theta - delta
    diff_dist = abs(diff_theta) * Ra
    p += diff_theta / ddelta
  end
  return tt, n, vert
end

# === Travel-time correction
function ttcorrection(pxpx,pxpy,pxpz,xducer_height,z,v,nz_st::Int64,numz::Int64,XDUCER_DEPTH,lat)
  tt_vert0, v_deep, v_range = ttvert0(pxpz,xducer_height,XDUCER_DEPTH,z,v,nz_st,numz)
  # --- Make synthetic data
  #xy_range = sqrt(pxpx^2 + pxpy^2)
  #xrange = xy_range*2.5
  #yrange = xy_range*2.5
  xrange = abs(pxpz)*2.5
  yrange = abs(pxpz)*2.5
  hrange = 5.0
  NN = 10000
  NP = 18
  xs = xrange*(rand(NN).-0.5)
  ys = yrange*(rand(NN).-0.5)
  hs = hrange*(rand(NN).-0.5) .+ xducer_height
  # --- Set inversion matrix
  d = zeros(NN)
  H = zeros(NN,NP)
  Rg, Rl = localradius(lat)
  for i in 1:NN
    tt_exact, Nitr, vert = xyz2tt(pxpx,pxpy,pxpz,xs[i],ys[i],hs[i],z,v,nz_st,numz,Rg,XDUCER_DEPTH)
    tt_sphere, theta = ttsphere(pxpx,pxpy,pxpz,xs[i],ys[i],hs[i],Rg,tt_vert0,v_deep,v_range)
    dh = hs[i] - xducer_height
    d[i] = tt_exact - tt_sphere
    H[i,1] = 1.0
    H[i,2] = theta
    H[i,3] = theta^2.0
    H[i,4] = theta^3.0
    H[i,5] = theta^4.0
    H[i,6] = theta^5.0
    H[i,7] = theta^6.0
    H[i,8] = theta^7.0
    H[i,9] = theta^8.0
    H[i,10] = dh
    H[i,11] = dh*theta
    H[i,12] = dh*theta^2.0
    H[i,13] = dh*theta^3.0
    H[i,14] = dh*theta^4.0
    H[i,15] = dh*theta^5.0
    H[i,16] = dh*theta^6.0
    H[i,17] = dh*theta^7.0
    H[i,18] = dh*theta^8.0
  end
  # --- Do inversion
  a = inv(transpose(H)*H)*transpose(H)*d
  dc = H*a
  dr = d - dc
  rms = std(dr)
  return tt_vert0, v_deep, v_range, a, rms
end

# === Calculate tt_appr
function xyz2tt_rapid(pxpx,pxpy,pxpz,rsx,rsy,rsz,Rg,tt_vert0,v_deep,v_range,xducer_height,cc0)
  delta_z = rsz - pxpz
  theta = sqrt((rsx-pxpx)^2.0+(rsy-pxpy)^2.0) / Rg
  b = Rg + rsz
  c = Rg + pxpz
  a = sqrt(b^2.0+c^2.0-2.0*b*c*cos(theta))
  vert = delta_z / a
  delta_h = delta_z - v_range
  tt_vert = tt_vert0 + delta_h / v_deep
  tt_sphere = tt_vert*a/delta_z
  tt_appr = tt_sphere
  dh = rsz - xducer_height
  for i in 1:9
    tt_appr += cc0[i]*theta^(i-1)
    tt_appr += cc0[i+9]*dh*theta^(i-1)
  end
  return tt_appr, tt_sphere, vert
end

# === Calculate tt_appr
function xyz2ttg_rapid(pxpx,pxpy,pxpz,rsx,rsy,rsz,Rg,tt_vert0,v_deep,v_range,xducer_height,cc0)
  delta_z = rsz - pxpz
  theta = sqrt((rsx-pxpx)^2.0+(rsy-pxpy)^2.0) / Rg
  xx = pxpx - rsx
  yy = pxpy - rsy
  xy_range = sqrt(xx^2 + yy^2)
  cote = xy_range / delta_z
  phi = atan(yy, xx)
  hh1 = cote*cos(phi)
  hh2 = cote*sin(phi)
  b = Rg + rsz
  c = Rg + pxpz
  a = sqrt(b^2.0+c^2.0-2.0*b*c*cos(theta))
  vert = delta_z / a
  delta_h = delta_z - v_range
  tt_vert = tt_vert0 + delta_h / v_deep
  tt_sphere = tt_vert*a/delta_z
  tt_appr = tt_sphere
  dh = rsz - xducer_height
  for i in 1:9
    tt_appr += cc0[i]*theta^(i-1)
    tt_appr += cc0[i+9]*dh*theta^(i-1)
  end
  return tt_appr, vert, hh1, hh2
end
