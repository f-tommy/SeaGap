export mat2tuple, tuple2mat, unixsort, unixsort2
"""
    mat2tuple(x)

Convert matrix `x` into tuple

# Example
    y = mat2tuple([1 2 3; 4 5 6])
""" 
function mat2tuple(x)
  num1 = size(x)[1]
  num2 = size(x)[2]
  y = []
  for n in 1:num1
    push!(y,(x[n,1:num2]))
  end
  return y
end

"""
    tuple2mat(y)

Convert tuple `y` into matrix

# Example
    x = tuple2mat((1,2))
"""
function tuple2mat(y)
  num1 = length(y)
  num2 = length(y[1])
  x = zeros(num1,num2)
  for n in 1:num1
    for m in 1:num2
      x[n,m] = (y[n])[m]
    end
  end
  return x
end

"""
   unixsort(x,k)

Sorting matrix `x` as with sort command of Unix system.
`k` is the keying colomun.

# Example
   x_new = unixsort(x,k)
"""
function unixsort(x,k::Int64)
  y = mat2tuple(x)
  sort!(y, by = col -> col[k])
  z = tuple2mat(y)
  return z
end

"""
   unixsort2(x,k1,k2)

Sorting matrix `x` as with sort command of Unix system.
`k1` is the primary keying colomun.
`k2` is the secoundary keying colomun.

# Example
   x_new = unixsort2(x,k1,k2)
"""
function unixsort2(x,k1::Int64,k2::Int64)
  y = mat2tuple(x)
  sort!(y, by = col -> (col[k1],col[k2]))
  z = tuple2mat(y)
  return z
end
