#using DelimitedFiles
#using Optim

export scale_est
"""
    scale_est(x0,s)

Scaling function for the long-term NTD parameters in `pos_array_mcmcpvg()` or `pos_array_mcmcpvgc()`.
`scale_est(x0,s)` transforms the original scale parameter `x0` into the scaled parameter `x`.

* `x0`: Original parameter (Vector)
* `s`: Scaling factor (Vector)

# Example
    x = scale_est([1.e-7],[-6.0])
"""
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
