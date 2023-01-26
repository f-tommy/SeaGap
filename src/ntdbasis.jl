# === Make temporal points for bspline interval
export mktbasis, tbspline3, retrieveb

"""
    mktbasis(NPB,t1,t2,NN)

Make temporal B-spline bases for modeling NTD.
`NPB` number of B-spline bases are deployed with a temporally-uniform interval.

Arguments:
* `NPB`: Number of B-spline bases
* `t1`: Time when an acoustic signal transmitting [sec] (vector)
* `t2`: Time when an acoustic signal recieving [sec] (vector)
* `NN`: Number of shots (norm of `t1` and `t2`)

Output:
* `smin`: Minimum values among `t1` and `t2`
* `smax`: Maximum values among `t1` and `t2`
* `ds`: Temporal interval of B-spline bases
* `tb`: Time at nodes that B-spine bases are deployed

# Example
    smin, smax, ds, tb = mktbasis(NPB,t1,t2,num)
"""
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
"""
    tbspline3(t,ds,tb,b,NPB)

Calculate value at `t` when providing coefficients `b` to temporal B-spline bases formed by `tb`, `ds`, and `NPB`.
                                                                                
Arguments:
* `t`: Time [sec]
* `ds`: Temporal interval of B-spline bases
* `tb`: Time at nodes that B-spine bases are deployed (vector)
* `b`:  Coeffcients for B-spline bases (vector)
* `NPB`: Number of B-spline bases

# Example
    tbs = tbspline3(t,ds,tb,b,NPB)
"""
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
"""
    retrieveb(NPB,tb,ds,t1,t2,num)

Identify unusable nodes of B-spline bases

Arguments:
* `NPB`: Number of B-spline bases
* `tb`: Time at nodes that B-spine bases are deployed (vector)
* `ds`: Temporal interval of B-spline bases

Output:
* `NPBV`: Number of usable B-spline bases
* `id`: if B-spline basis is not used, `id=0` (vector with norm of `NPB`)

# Example
    NPBV, id = retrieveb(NPB,tb,ds,t1,t2,num) 
"""
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
