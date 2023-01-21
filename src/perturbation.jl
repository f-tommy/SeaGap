#using Random

export perturbation_single
function perturbation_single(x0,dx,x1,x2)
  dr = (rand()-0.5)*2
  x = x0 + dr*dx
  if x < x1 || x > x2
    x = x0
  end
  return x
end
