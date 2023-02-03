#using LinearAlgebra
#using Dierckx

export interpolate_gps
"""
    interpolate_gps(t0,x0,y0,z0,h0,p0,r0,tp0,tt0,nump; ks)

Iterpolate GNSS position and attitude when an acoustic signal is transmitted from a sea-surface transducer and arrived at the sea-surface transducer

* `t0`: Time [sec] (`nump` size of vector)
* `x0`: EW position [m] (`nump` size of vector)
* `y0`: NS position [m] (`nump` size of vector)
* `z0`: UD position [m] (`nump` size of vector)
* `h0`: Heading [deg] (`nump` size of vector)
* `p0`: Pitch [deg] (`nump` size of vector)
* `r0`: Roll [deg] (`nump` size of vector)
* `tp0`: Time when an coustic signal transmitted [sec] (`numk` size of vector)
* `tt0`: Travel-time [sec] (`tp0`+`tt0`->an acoustic signal is arrived: `numk` size of vector)
* `nump`: Number of data for `t0`...`r0`
* `numk`: Number of transponders for `tp0` and `tt0`
* `ks`: Degree of spline function for interpoaltion (`ks=3` by default)

# Example
    xd0,yd0,zd0,hd0,pd0,rd0,xd1,yd1,zd1,hd1,pd1,rd1,tp1 = interpolate_gps(t0,x0,y0,z0,h0,p0,r0,tp0,tt0,nump)
"""
function interpolate_gps(t0,x0,y0,z0,h0,p0,r0,tp0,tt0,nump,numk::Int64,maxp::Int64; ks=3)
  # --- Make cubic spline function
  cubic_x = Spline1D(t0, x0, k=ks)
  cubic_y = Spline1D(t0, y0, k=ks)
  cubic_z = Spline1D(t0, z0, k=ks)
  cubic_h = Spline1D(t0, h0, k=ks)
  cubic_p = Spline1D(t0, p0, k=ks)
  cubic_r = Spline1D(t0, r0, k=ks)
  # --- initialize
  tp1 = zeros(numk,maxp)
  xd0 = zeros(numk,maxp)
  yd0 = zeros(numk,maxp)
  zd0 = zeros(numk,maxp)
  xd1 = zeros(numk,maxp)
  yd1 = zeros(numk,maxp)
  zd1 = zeros(numk,maxp)
  hd0 = zeros(numk,maxp)
  pd0 = zeros(numk,maxp)
  rd0 = zeros(numk,maxp)
  hd1 = zeros(numk,maxp)
  pd1 = zeros(numk,maxp)
  rd1 = zeros(numk,maxp)
  for k in 1:numk
    kk = Int(nump[k])
    tp1[k,1:kk] = tp0[k,1:kk] + tt0[k,1:kk]
    xd0[k,1:kk] = cubic_x.(tp0[k,1:kk])
    yd0[k,1:kk] = cubic_y.(tp0[k,1:kk])
    zd0[k,1:kk] = cubic_z.(tp0[k,1:kk])
    xd1[k,1:kk] = cubic_x.(tp1[k,1:kk])
    yd1[k,1:kk] = cubic_y.(tp1[k,1:kk])
    zd1[k,1:kk] = cubic_z.(tp1[k,1:kk])
    hd0[k,1:kk] = cubic_h.(tp0[k,1:kk])
    pd0[k,1:kk] = cubic_p.(tp0[k,1:kk])
    rd0[k,1:kk] = cubic_r.(tp0[k,1:kk])
    hd1[k,1:kk] = cubic_h.(tp1[k,1:kk])
    pd1[k,1:kk] = cubic_p.(tp1[k,1:kk])
    rd1[k,1:kk] = cubic_r.(tp1[k,1:kk])
  end
  return xd0,yd0,zd0,hd0,pd0,rd0,xd1,yd1,zd1,hd1,pd1,rd1,tp1
end

