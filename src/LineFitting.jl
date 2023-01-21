#using Statistics
#using LinearAlgebra
#using Distributions

struct LineFitting{T<:Union{Matrix{<:Real},Vector{<:Real}}}
  coef::T
  coefstd::T
  pred::T
  misfit::T
  sigma2::Float64
  pred_new::T
  pred_lower::T
  pred_upper::T
  alpha::Float64
  newX::T
end

"""
   linefit(x,y,w; alpha,newX)

Perform line fitting: `y` = A `x` + B
* `x`: Vector{Float64} 
* `y`: Vector{Float64} 
* `w`: Vector{Float64}, optional wieghts of `y`
* alpha: Float64, Percentage of Condidence Interval
* `newX`: Vector{Float64}, calculate predicted values for `newX`

Output:
* lf.coef: estimates of A and B
* lf.coefstd: Standard deviations of A and B
* lf.pred: Predicted values for `x`
* lf.misfit: Misfits
* lf.sigma2: Std. of misfits
* lf.pred_new: Predicted values for `newX` if `newX` is given
* lf.pred_lower: Lower `alpha`% CI of predictions for `newX`
* lf.pred_upper: Upper `alpha`% CI of predictions for `newX`
* lf.alpha: Return `alpha`
* lf.newX: Return `newX`

# Example

lf = linefit(x,y)

"""

export linefit
function linefit(x,y,w0=[];alpha=0.95::Float64,newX=[])
  # --- Check the input values
  nx = length(x); ny = length(y); nw = length(w0); nnx = length(newX)
  if nx != ny
    error("Set vectors with the same length for x and y")
  end
  if nx < 3
    error("Number of observations must be >= 3")
  end
  if nw > 0 && nw != nx
    error("Length of vector w must be same with those of x and y")
  end
  # --- newX
  if nnx == 0
    newX = copy(x)
  end
  # --- Equal weights
  if nw == 0
    w0 = ones(nx)
  end
  # --- Estimation
  d = copy(y); H = ones(nx,2); W = diagm(w0)
  for n in 1:nx
    H[n,2] = x[n]
  end
  # --- Inversion
  HWd = transpose(H)*W*d; HWH = transpose(H)*W*H
  Hinv = inv(HWH)
  a = Hinv*HWd
  # --- calculation and residuals
  pred = H*a
  misfit = d - pred
  # --- Obtain covariance matrix
  sa = dot(misfit,W*misfit)
  sigma2 = sa/(nx-2)
  C = sigma2*Hinv
  coefstd = sqrt.(diag(C))
  coef = a
  # --- Confidence interval
  nnx = length(newX)
  Hn = ones(nnx,2)
  for n in 1:nnx
    Hn[n,2] = newX[n]
  end
  pred_new = Hn*a
  ta = quantile.(Distributions.TDist(nx-2),[(1-alpha)/2,1-(1-alpha)/2])
  ones(nnx)
  Hn*Hinv*transpose(Hn)
  pred_lower = pred_new + ta[1]*sqrt.(diag(Hn*C*transpose(Hn)))
  pred_upper = pred_new + ta[2]*sqrt.(diag(Hn*C*transpose(Hn)))
  return LineFitting(coef,coefstd,pred,misfit,sigma2,pred_new,pred_lower,pred_upper,alpha,newX)
end
