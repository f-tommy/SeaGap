#using Random

export perturbation_single
"""
    perturbation_single(x0,dx,x1,x2)

Add perturbation to `x0`.
The pertubation is obtained from uniform distribution with a width of `dx`.
If new value is out of (`x1`,`x2`), the original value `x0` is returned.

# Example
    perturbation_single(5.1,0.1,0.0,20.0)
"""
function perturbation_single(x0,dx,x1,x2; dr=(rand()-0.5)*2)
  x = x0 + dr*dx
  if x < x1 || x > x2
    x = x0
  end
  return x
end
