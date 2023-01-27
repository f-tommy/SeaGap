export scale
"""
    scale(x,s)

Scaling function for the long-term NTD parameters in `pos_array_mcmcpvg()` or `pos_array_mcmcpvgc()`.
`scale(x,s)` transform the sclaed paramter `x` into the original scale parameter `y`.

* `x`: Scaled parameter
* `s`: Scaling factor

# Example
    y = scale(0.1,-7.0)
"""
function scale(x,s)
  y = x*10.0^(abs(x)+s)
  return y
end
