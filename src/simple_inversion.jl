#using Dates
#using LinearAlgebra

export simple_inversion
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
