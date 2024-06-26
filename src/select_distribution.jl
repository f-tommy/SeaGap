export select_dist_gmode, select_dist_dmode, select_dist_rmode
export select_dist_hmode

function select_dist_gmode(gmode,gm,gs)
  if gmode == 0
    gdist = Distributions.Uniform(gm-gs,gm+gs)
  elseif gmode == 1
    gdist = Distributions.Normal(gm,gs)
  elseif gmode == 2
    gdist = Distributions.Cauchy(gm,gs)
  elseif gmode == 3
    gdist = Distributions.Laplace(gm,gs)
  end
  return gdist
end

function select_dist_dmode(dmode,dm,ds,dep)
  if dmode == 0
    ddist = Distributions.Uniform(dm-ds,dm+ds)
  elseif dmode == 1
    ddist = Distributions.Uniform(0.0,dep)
  elseif dmode == 2
    ddist = Distributions.Normal(dm,ds)
  elseif dmode == 3
    ddist0 = Distributions.Normal(dm,ds)
    ddist = Distributions.truncated(ddist0,0.0,dep)
  elseif dmode == 4
    ddist = Distributions.Cauchy(dm,ds)
  elseif dmode == 5
    ddist0 = Distributions.Cauchy(dm,ds)
    ddist = Distributions.truncated(ddist0,0.0,dep)
  elseif dmode == 6
    ddist = Distributions.Laplace(dm,ds)
  elseif dmode == 7
    ddist0 = Distributions.Laplace(dm,ds)
    ddist = Distributions.truncated(ddist0,0.0,dep)
  end
  return ddist
end

function select_dist_rmode(rmode,rm,rs)
  if rmode == 0
    rdist = Distributions.Uniform(rm-rs,rm+rs)
  elseif rmode == 1
    rdist = Distributions.Normal(rm,rs)
  elseif rmode == 2
    rdist = Distributions.Cauchy(rm,rs)
  elseif rmode == 3
    rdist = Distributions.Laplace(rm,rs)
  end
  return rdist
end

function select_dist_hmode(hmode,hs)
  if hmode == 0
    hdist = Distributions.Uniform(0.0,100.0)
  elseif hmode == 1
    hdist0 = Distributions.Normal(0.0,10^hs)
    hdist = Distributions.truncated(hdist0; lower=0.0)
  elseif hmode == 2
    hdist0 = Distributions.Cauchy(0.0,10^hs)
    hdist = Distributions.truncated(hdist0; lower=0.0)
  end
  return hdist
end
