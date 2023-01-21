#using Statistics
export runmed, runave
function runmed(x,k::Int64)
  if k < 3
    error("runmed: k must be >=3")
  end
  if k%2 == 0
    error("runmed: k must be odd")
  end
  num=length(x)
  if num < k
    error("runmed: Length of input data vector should be >= k")
  end

  kk = Int(floor(k/2))
  numk = num - kk
  y = copy(x)
  for i in 1:num
    d = zeros(k)
    if i <= kk
      d[1:k] = x[1:k]
    elseif i > numk
      d[1:k] = x[num-k+1:num]
    else
      d[1:k] = x[i-kk:i+kk]
    end
    y[i] = median(d)
  end
  return y
end

function runave(x,k::Int64)
  if k < 3
    error("runave: k must be >=3")
  end
  if k%2 == 0
    error("runave: k must be odd")
  end
  num=length(x)
  if num < k
    error("runmed: Length of input data vector should be >= k")
  end

  kk = Int(floor(k/2))
  numk = num - kk
  y = copy(x)
  for i in 1:num
    d = zeros(k)
    if i <= kk
      d[1:k] = x[1:k]
    elseif i > numk
      d[1:k] = x[num-k+1:num]
    else
      d[1:k] = x[i-kk:i+kk]
    end
    y[i] = mean(d)
  end
  return y
end
