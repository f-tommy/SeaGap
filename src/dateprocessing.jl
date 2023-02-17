#using Dates
#using Printf

export date2cal, date2cal_txt, date2year, date2year_txt, date2sec, date2sec_txt, sec2date, sec2date_txt, sec2year, sec2year_txt

"""
    date2cal(t)

String type date `t` is converted into cumulative days `cal` in `year`

# Example
    year, cal = date2cal("2020-10-10")
"""
function date2cal(t)
  tv = Dates.DateTime(t) # Convert into Date-type
  year = Dates.year(tv)
  cal = Dates.dayofyear(tv)
  return year,cal
end

"""
    date2cal_txt(fn0,fn)                                                                 

Read text file `fn0`(file name), convert date at 1st column into cumulative days, and write `fn` as an output file.

# Example
    date2cal_txt(fn0,fn)
"""

function date2cal_txt(fn0::String,fn::String)
  f = open(fn0)
  a = readlines(f)
  close(f)
  num1 = size(a)[1]
  num2 = size(split(a[1]))[1]
  open(fn,"w") do out
  if num2 >= 2
    for i in 1:num1
      y, c = date2cal(split(a[i])[1])
      d = join(split(a[i])[2:num2]," ")
      @printf(out,"%4d %03d %s\n",y,c,d)
    end
  else
    for i in 1:num1
      y, c = date2cal(a[i])
      @printf(out,"%4d %03d\n",y,c)
    end
  end
  end
end

"""
    date2year(t)

String type date `t` is converted into `year`

# Example
    year = date2year("2020-10-10T00:10:21")
"""
function date2year(t)
  tv = Dates.DateTime(t) # Convert into Date-type
  year0 = Dates.year(tv)
  cal = Dates.dayofyear(tv)
  aday = Dates.daysinyear(tv)
  year = year0 + cal/aday
  return year
end

"""
    date2year_txt(fn0,fn)

Read text file `fn0`(file name), convert date at 1st column into year, and write `fn` as an output file.

# Example
    date2year_txt(fn0,fn)
"""
function date2year_txt(fn0::String,fn::String)
  f = open(fn0)
  a = readlines(f)
  close(f)
  num1 = size(a)[1]
  num2 = size(split(a[1]))[1]
  open(fn,"w") do out
  if num2 >= 2
    for i in 1:num1
      y = date2year(split(a[i])[1])
      d = join(split(a[i])[2:num2]," ")
      @printf(out,"%4.12f %s\n",y,d)
    end
  else
    for i in 1:num1
      y = date2year(a[i])
      @printf(out,"%4.12f\n",y,c)
    end
  end
  end
end

"""
    date2sec(t,t0)

String type date `t` is converted into cumulative seconds `sec` from `t0`.
`t0` is set to "2000-01-01T12:00:00" by default.

# Example
    sec = date2sec("2020-10-10T00:10:21")
"""
function date2sec(t,t0="2000-01-01T12:00:00")
  t0 = Dates.DateTime(t0)
  tv0 = split(t,'.')
  numk = size(tv0)[1]
  if numk == 1
    tv = Dates.datetime2unix(Dates.DateTime(tv0[1]))
  elseif numk == 2
    tv = Dates.datetime2unix(Dates.DateTime(tv0[1])) + parse(Float64,string("0.",tv0[2]))
  else
    error("date2sec.jl: Format is collapsed")
  end
  return tv - Dates.datetime2unix(t0)
end

"""
    date2sec_txt(fn0,fn,t0)                                                                 

Read text file `fn0`(file name), convert date at 1st column into cumulative seconds from `t0`, and write `fn` as an output file.
`t0` is set to "2000-01-01T12:00:00" by default.

# Example
    date2sec_txt(fn0,fn)
"""
function date2sec_txt(fn0::String,fn::String,t0="2000-01-01T12:00:00")
  t0 = Dates.DateTime(t0)
  f = open(fn0)
  a = readlines(f)
  close(f)
  num1 = size(a)[1]
  num2 = size(split(a[1]))[1]
  open(fn,"w") do out
  if num2 >= 2
    for i in 1:num1
      dt = date2sec(split(a[i])[1],t0)
      d = join(split(a[i])[2:num2]," ")
      @printf(out,"%10.6f %s\n",dt,d)
    end
  else
    for i in 1:num1
      dt = date2sec(a[i],t0)
      @printf(out,"%10.6f\n",dt)
    end
  end
  end
end

"""
    sec2date(j,t0; digits)

Cumulative seconds `j` from `t0` is converted into `date`.
`t0` is set to "2000-01-01T12:00:00" by default.
`digits` is the digits of seconds; 6 by default.

# Example
   date  = sec2date(12345.678)
"""
function sec2date(j,t0="2000-01-01T12:00:00"; digits=6::Int64)
  t0 = Dates.DateTime(t0)
  j0 = Dates.datetime2unix(t0)
  jj = j + j0
  jj1 = floor(jj)
  jj2 = split(string(round(jj-floor(jj),digits=digits)),'.')
  dat0 = string(Dates.unix2datetime(jj1))
  dat = string(dat0,".",jj2[2])
  return dat 
end

"""
    sec2date_txt(fn0,fn,t0; digits)

Read text file `fn0`(file name), convert cumulative seconds at 1st column into date, and write `fn` as an output file.
`t0` is set to "2000-01-01T12:00:00" by default.
`digits` is the digits of seconds; 6 by default.

# Example
    sec2date_txt(fn0,fn)
"""
function sec2date_txt(fn0,fn,t0="2000-01-01T12:00:00";digits=6::Int64)
  f = open(fn0)
  a = readlines(f)
  close(f)
  num1 = size(a)[1]
  num2 = size(split(a[1]))[1]
  open(fn,"w") do out
  if num2 >= 2
    for i in 1:num1
      j = sec2date(parse(Float64,split(a[i])[1]),t0,digits=digits)
      d = join(split(a[i])[2:num2]," ")
      @printf(out,"%s %s\n",string(j),d)
    end
  else
    for i in 1:num1
      j = sec2date(parse(Float64,a[i]),t0,digits=digits)
      @printf(out,"%s\n",string(j))
    end
  end
  end
end

"""
    sec2year(j,t0)

Cumulative seconds `j` from `t0` is converted into `year`.
`t0` is set to "2000-01-01T12:00:00" by default.

# Example
   year  = sec2year(12345.678)
"""
function sec2year(j,t0="2000-01-01T12:00:00")
  dat0 = sec2date(j,t0)
  t = split(dat0,'.')
  numk = size(t)[1]
  if numk == 1
    dat = dat0
  else
    dat = t[1]
  end
  year = date2year(dat)
  return year
end

"""
    sec2year_txt(fn0,fn,t0)

Read text file `fn0`(file name), convert cumulative seconds at 1st column into year, and write `fn` as an output file.
`t0` is set to "2000-01-01T12:00:00" by default.

# Example
    sec2year_txt(fn0,fn)
"""
function sec2year_txt(fn0::String,fn::String,t0="2000-01-01T12:00:00")
  f = open(fn0)
  a = readlines(f)
  close(f)
  num1 = size(a)[1]
  num2 = size(split(a[1]))[1]
  open(fn,"w") do out
  if num2 >= 2
    for i in 1:num1
      j = sec2year(parse(Float64,split(a[i])[1]),t0)
      d = join(split(a[i])[2:num2]," ")
      @printf(out,"%4.12f %s\n",j,d)
    end
  else
    for i in 1:num1
      j = sec2year(parse(Float64,a[i]),t0)
      @printf(out,"%4.12f\n",j)
    end
  end
  end
end
