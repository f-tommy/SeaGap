#=
  Functions for anttena2tr
    Euler Rotation 
  written by Fumiaki Tomita 2022/06/29
=#
#using LinearAlgebra
export rotate, anttena2tr

"""
    rotate(x,ex,ey,ez)

Euler roation for three component vector `x` depending on `ex`, `ey`, and `ez`.

# Example
    y  = rotate(x,ex,ey,ez)
"""
function rotate(x,ex,ey,ez)
  a = sqrt(ex^2+ey^2+ez^2)
  if abs(a) < 1.e-10
    y = x
  else
    r = zeros(3,3)
    cosa = cos(a)
    sina = sin(a)
    r[1,1] = ex^2 /a^2*(1-cosa) + cosa
    r[1,2] = ex*ey/a^2*(1-cosa) - ez/a*sina
    r[1,3] = ex*ez/a^2*(1-cosa) + ey/a*sina
    r[2,1] = ex*ey/a^2*(1-cosa) + ez/a*sina
    r[2,2] = ey^2 /a^2*(1-cosa) + cosa
    r[2,3] = ey*ez/a^2*(1-cosa) - ex/a*sina
    r[3,1] = ex*ez/a^2*(1-cosa) - ey/a*sina
    r[3,2] = ey*ez/a^2*(1-cosa) + ex/a*sina
    r[3,3] = ez^2 /a^2*(1-cosa) + cosa
    y = r * x
  end
  return y
end

"""
    anttena2tr(x0,y0,z0,h,p,r,e)
                                                                                  
Transform GNSS antenna position in an locally orthogonal coordinate (x0,y0,z0) into sea-surface transducer position (`x`,`y`,`z`) depending on vessel's attitude (heading:`h`, pitch:`p`, roll:`r`) with three component of antenna-transducer offsets `e` (the transduecer position from the antenna position).

# Example
    x, y, z = anttena2tr(x0,y0,z0,h,p,r,e)
"""
function anttena2tr(x0,y0,z0,h,p,r,e)
  # --- deg2rad
  hh = h*pi/180.0
  pp = p*pi/180.0
  rr = r*pi/180.0
  # --- heading transform
  a = -hh
  ex = 0.0
  ey = 0.0
  ez = a
  f = rotate(e,ex,ey,ez)
  # --- pitch transform
  a = pp
  ex = a*sin(pi/2.0)*cos(2*pi-hh)
  ey = a*sin(pi/2.0)*sin(2*pi-hh)
  ez = a*cos(pi/2.0)
  g = rotate(f,ex,ey,ez)
  # --- roll transform
  a = rr
  ex = a*sin(pi/2.0-pp)*cos(pi/2.0-hh)
  ey = a*sin(pi/2.0-pp)*sin(pi/2.0-hh)
  ez = a*cos(pi/2.0-pp)
  f = rotate(g,ex,ey,ez)
  # --- transform
  x = x0 + f[1]
  y = y0 + f[2]
  z = z0 + f[3]
  return x, y, z
end
