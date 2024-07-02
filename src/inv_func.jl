export smoothing_matrix_1d
function smoothing_matrix_1d(NP)
  if NP >= 2
    Lg = zeros(NP,NP)
    Lg[1,1] = 1.0; Lg[1,2] = -1.0
    Lg[NP,NP] = 1.0; Lg[NP,NP-1] = -1.0
    if NP >= 3
      for i in 2:NP-1
        Lg[i,i] = 2.0
        Lg[i,i-1] = -1.0
        Lg[i,i+1] = -1.0
      end
    end
  else
    Lg = ones(1)
  end
  return Lg
end
