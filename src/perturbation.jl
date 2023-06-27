#using Random

export perturbation_param
export perturbation_param_nolimit
export perturbation_along
"""
    perturbation_param(x0,dx,x1,x2)

Add perturbation to `x0`.
The pertubation is obtained from uniform distribution with a width of `dx`.
If new value is out of (`x1`,`x2`), the original value `x0` is returned.

# Example
    perturbation_param(5.1,0.1,0.0,20.0)
"""
function perturbation_param(x0,dx,x1,x2; dr=(rand()-0.5)*2)
  x = x0 + dr*dx
  if x < x1 || x > x2
    x = x0
  end
  return x
end

function perturbation_param_nolimit(x0,dx; dr=(rand()-0.5)*2)
  x = x0 + dr*dx
  return x
end

function perturbation_along(x0,dx,x1,x2,y0,dy,y1,y2)
  dr = (dx + dy)/2.0
  delta = (rand()-0.5)*2.0*dr
  x = x0 + delta*rand()
  y = y0 + delta*rand()
  if x < x1 || x > x2
    x = x0
  end
  if y < y1 || y > y2
    y = y0
  end
  return x,y
end
