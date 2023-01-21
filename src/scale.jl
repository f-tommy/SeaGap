export scale
function scale(x,s)
  y = x*10.0^(abs(x)+s)
  return y
end
