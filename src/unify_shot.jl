export unify_shot
function unify_shot(tp0,nump,numk::Int64)
  t_all = []
  for k in 1:numk
    append!(t_all,tp0[k,1:Int(nump[k])])
  end
  tp = sort(unique(t_all))
  num = size(tp)[1]
  println(stderr,"     Unified shot number: $num")
  id = zeros(num); kd = zeros(numk,num) 
  for n in 1:num
    for k in 1:numk
      id[n] += sum(tp[n] .== tp0[k,1:Int(nump[k])])
      for i in 1:Int(nump[k])
        if tp0[k,i] == tp[n]
          kd[k,n] = i
        end
      end
    end
  end
  return num,Int.(tp),Int.(id),Int.(kd)
end

