using Memoization
using Nemo



# add "a" to row k of tab, return nothing if fails
function add_tabi(tab, k, a)
    len = length(tab)

    if k > len + 1
        return nothing
    end

    if k == 1
        if a >= last(tab[1])
            temp = [[tab[1]; a], tab[2:end]...]
            return temp
        else
            return nothing
        end
    end

    if k > 1 && k < len + 1
        lam = length(tab[k])
        if lam < length(tab[k-1]) && 
           a >= last(tab[k]) && 
           a > tab[k-1][lam+1]
            temp = [tab[1:k-1]..., [tab[k]; a], tab[k+1:end]...]
            return temp
        else
            return nothing
        end
    end

    if k == len + 1
        if a > first(tab[len])
            temp = [tab..., [a]]
            return temp
        else
            return nothing
        end
    end
end



@memoize function ssyt(la, ns; ff = fill(maximum(ns), length(la)), mu = fill(0, length(la)))
    len = length(la)
    mmu = copy(mu)
    
    while length(mmu) < len
        push!(mmu, 0)
    end
    
    if length(mmu) > len || mmu[end] > la[end]
        println("error: mu not contained in lambda")
        return nothing
    end
    
    tabs = []
    
    if len == 0
        return []
    end
    
    if len == 1
        if la[1] == mmu[1]
            push!(tabs, [fill(0, mmu[1])] )
        end
        if la[1] == mmu[1] + 1
            for i in 1:length(ns)
                if ns[i] <= ff[1]
                    push!(tabs, [append!( copy(fill(0, mmu[1])), ns[i])] )
                end
            end
            return tabs
        end
        if la[1] > mmu[1] + 1
            tabs1 = ssyt([la[1] - 1], ns, ff=ff, mu=mmu)
            for tt in tabs1
                for k in 1:length(ns)
                    if ns[k] <= ff[1]
                        t2=add_tabi(tt,1,ns[k])
                        if t2!=nothing
                          push!(tabs, t2)
                        end
                    end
                end
            end
        end
        return tabs
    end
    
    if len > 1
        la1 = la[1:len-1]
        mu1 = mmu[1:len-1]
        
        tabs1 = ssyt(la1, ns, ff=ff, mu=mu1)
        
        if mmu[end] > 0
            tabs2 = []
            for tt in tabs1
                push!(tabs2, vcat( tt, [fill(0, mmu[len])] ) )
            end
            tabs1 = tabs2
        end
        
        for i in 1:(la[len] - mmu[len])
            tabs2 = []
            for tt in tabs1
                for k in 1:length(ns)
                    if ns[k] <= ff[len]
                        t2 = add_tabi(tt, len, ns[k])
                        if t2!=nothing
                          push!(tabs2, t2)
                        end
                    end
                end
            end
            tabs1 = tabs2
        end
        tabs = tabs1
    end
        
    return tabs
end




######
struct DoublePolyRing
           ring::ZZMPolyRing
           x_vars::Vector{ZZMPolyRingElem}
           y_vars::Vector{ZZMPolyRingElem}
       end

#####
function xy_ring(n,m)
# construct polynomial ring in n xvars and m yvars
  local xvars = ["x$(i)" for i=1:n]
  local yvars = ["y$(i)" for i=1:m]

  R,all_vars = polynomial_ring(ZZ,vcat(xvars,yvars))

  x = all_vars[1:n]
  y = all_vars[n+1:n+m]

  return DoublePolyRing(R,x,y), x, y

end




function tab2bin( tab, RR::DoublePolyRing; xoffset = 0, yoffset = 0 )
# return the product of binomials
  len = length(tab)

  x = RR.x_vars
  y = RR.y_vars

  n = length(x)
  m = length(y)

  bin = RR.ring(1)

  for i=1:len
    for j=1:length( tab[i] )
      tt = tab[i][j]
      if tt>0
        p=RR.ring(0)
        if tt+xoffset<=n
          p=p+x[tt+xoffset]
        end
        if tt+j-i+yoffset<=m && tt+j-i+yoffset>0
          p=p+y[tt+j-i+yoffset]
        end
      end
      bin = bin*p
    end
  end

  return bin

end


function ssyt2pol( tabs, RR::DoublePolyRing; xoffset=0, yoffset=0 )
# sum of binomials for a set of tableaux

  pol=RR.ring(0)

  for tab in tabs
    pol = pol + tab2bin( tab, RR; xoffset=xoffset, yoffset=yoffset )
  end

  return pol

end


function schur_poly( la, ff::Vector{Int}, RR::DoublePolyRing=xy_ring( length(la) , length(la)+la[1] )[1]; mu = fill(0, length(la)), xoffset=0, yoffset=0 )
  if length(la)==0
    return RR.ring(1)
  end

  tbs = ssyt( la, collect(1:maximum(ff)), ff=ff, mu=mu )

  pol = ssyt2pol( tbs, RR; xoffset=xoffset, yoffset=yoffset )

end



function schur_poly( la, ff::Int, RR::DoublePolyRing=xy_ring( length(la) , length(la)+la[1] )[1]; mu = fill(0, length(la)), xoffset=0, yoffset=0 )
  if length(la)==0
    return RR.ring(1)
  end

  tbs = ssyt( la, collect(1:ff), ff=fill(ff,ff), mu=mu )

  pol = ssyt2pol( tbs, RR; xoffset=xoffset, yoffset=yoffset )

end

