#using DelimitedFiles
#using Optim

export scale_est
function scale_est(d,a)
  num = size(d)[1]
  # --- Inversion
  y = zeros(num)
  for n in 1:num
    g(x) = sqrt(abs(scale(x[1],a[n]) - d[n]))
    x0 = [0.0]
    res = optimize(g,x0)
    y[n] = Optim.minimizer(res)[1]
  end
  return y
end
