# === Make temporal points for bspline interval
export mktbasis, tbspline3, retrieveb
function mktbasis(NPB::Int64,t1,t2,NN::Int64)
  if NPB < 1
    error("mktbasis: number of B-spline basis must be integer and over 1")
  end
  if NN < 1
    error("mktbasis: number of data must be over 1")
  end
  smin = findmin(t1)[1]
  smax = findmax(t2)[1]
  ds = (smax - smin) / (NPB - 3.0)
  println(stderr,"     B-spline interval:",ds)
  tb = zeros(NPB)
  for n in 1:NPB
    tb[n] = smin + ds*(n - 2.0)
  end
  println(stderr,"     Min sec:",smin," ",tb[1])
  println(stderr,"     Max sec:",smax," ",tb[NPB])
  return smin, smax, ds, tb
end 

# === Temporal B-spline function
function tbspline3(t,ds,tb,b,NPB::Int64)
  if NPB < 1
    error("tbspline3: number of B-spline basis must be integer and over 1")
  end
  a = 0.0
  for n in 1:NPB
    tb0 = sqrt(((t-tb[n])/ds)^2.0)
    x = abs(tb0)

    if x <= 1.0
      y = (3.0*x^3.0 - 6.0*x^2.0 + 4.0) / 6.0
    elseif x > 1.0 && x <= 2.0
      y = (2.0 - x)^3.0 / 6.0
    elseif x > 2.0
      y = 0.0
    end

    a += y * b[n]
  end
  return a
end

# === Retrieve B-spline knot
function retrieveb(NPB::Int64,tb,ds,t1,t2,num::Int64)
  if NPB < 1
    error("retrieveb: number of B-spline basis must be integer and over 1")
  end
  if num < 1
    error("retrieveb: number of data must be over 1")
  end
  NPBV = 0
  id = zeros(NPB)
  for m in 1:NPB
    b = zeros(NPB)
    b[m] = 1.0
    tmp1 = 0.0
    for n in 1:num
      tmp1 += tbspline3((t1[n]+t2[n])/2.0,ds,tb,b,NPB)
    end
    if tmp1 != 0.0
      NPBV += 1
      id[m] = NPBV
    else
      id[m] = 0
    end
  end
  println(stderr,"     B-spline knot: $NPB -> $NPBV")
  return Int(NPBV), Int.(round.(id))
end
