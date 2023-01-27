#using Dates
#using LinearAlgebra

export simple_inversion
"""
    simple_inversion(d,H)

Simple inversion by linear least squares method without any constraints.

* `d`: Data vector
* `H`: Karnel matrix
* `dcal`: Calculation vector
* `dres`: Residual vector (`d`-`dcal`)
* `a`: Estimated parameters
* `e`: iStd of the estimated parameters

# Example
   dcal, dres, a, e = simple_inversion(d,H)
"""
function simple_inversion(d,H)
  # --- Set number
  NN0 = length(d)
  NN = size(H)[1]
  NP = size(H)[2]
  if NN != NN0
    error(" simple_inversion: Dimension of d and row of H must be same")
  end
  if NN < NP
    error(" simple_inversion: number of rows must be over number of columns for H")
  end
  println(" Inversion matrix: ",NN,"*",NP)
  # --- Inversion
  Hd = transpose(H)*d
  HH = transpose(H)*H
  Hinv = inv(HH)
  a = Hinv*Hd
  # --- calculation and residuals
  dcal = H*a
  dres = d-dcal
  # --- Obtain covariance matrix
  sa = dot(dres,dres)
  sigma2 = sa/(NN-NP)
  C = sigma2*Hinv
  e = sqrt.(diag(C))
  # --- Output
  return dcal, dres, a, e
end
