# Tools for computing Schubert polynomials in Julia, using BPDs
# David Anderson, October 2023.


#Using Nemo

#include("ssyt.jl")
include("SSYT.jl")
include("BpdBase.jl")
include("Drifts.jl")

#=
#These are currently loaded in ssyt.jl
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

=#
######


function bpd2bin( bpd, R::DoublePolyRing=xy_ring( size(bpd.m)[1]-1, size(bpd.m)[2]-1 )[1]  )
# product of binomials for bpd
# requires DoublePolyRing
# can get single polyn by using no y_vars
  local n=size(bpd.m)[1]-1
  bin = R.ring(1)

  x = R.x_vars
  y = R.y_vars

  local aa=length(x)
  local bb=length(y)

  for i=1:n
    for j=1:n
      if bpd.m[i,j]==0
        p=R.ring(0)
        if i<=aa
          p=p+x[i]
        end
        if j<=bb
          p=p+y[j]
        end        
        bin = bin*p
      end
    end
  end

  return bin

end


function bpds2pol( bpds, R::DoublePolyRing )
# form polynomial from collection of bpds

  pol=R.ring(0)

  for bp in bpds
    pol = pol+bpd2bin(bp,R)
  end

  return(pol)

end


function schub_bpd( w, R::DoublePolyRing=xy_ring( length(w)-1, length(w)-1 )[1]  )
# compute schubert pol by bpd formula
  bpds=all_bpds(w)

  pol=R.ring(0)

  for bp in bpds
    pol = pol+bpd2bin(bp,R)
  end

  return(pol)

end




function schub_flats( w, R::DoublePolyRing=xy_ring( length(w)-1, length(w)-1 )[1] )
# compute schubert pol by drift class formula

  fbpds = flat_bpds(w)

  pol = R.ring(0)

  for b in fbpds
    b=markconfig(b)
    pol = pol+dc2sd(b,R)
  end

  return pol

end


function dc2sd( dc::Drift, R::DoublePolyRing=xy_ring( size(dc.m)[1]-1, size(dc.m)[2]-1 )[1]  )
# drift config to s-polynomial
# must take marked config as input

  local n=size(dc.m)[1]

  for k=2*n:-1:2
    for i=maximum([1,k-n]):minimum([n,k-1])
      if isa( dc.m[i,k-i], Tuple ) && dc.m[i,k-i][2]
        (dc1,dc2)=drift_split( dc, i, k-i )
        return ( dc2sd( dc1, R ) + dc2sd( dc2, R ) )
      end
    end
  end

  sd = R.ring(1)

  tcomps = tableau_components(dc)

  for tt in tcomps
    sd = sd*schur_poly( tt[1], tt[2], R; mu=tt[3], xoffset=tt[4][1], yoffset=tt[4][2], rowmin=true )
  end

  return sd

end


function tableau_components(dc)
# return labelled tableaux for a drift config dc

  local n=size(dc.m)[1]

  if !isflat(dc)
    return( tableau_components( nw_reset(dc) ) )
  end

  local lyds=Vector{Vector}([])

  local corners=Vector{Tuple{Int,Int}}([])

  for i=1:n
    for j=1:n
        if !( (i,j) in corners) && isa(dc.m[i,j],Tuple) && ((i,j)==(1,1) || (i>1 && j>1 && !isa( dc.m[i-1,j], Tuple) && !isa( dc.m[i,j-1],Tuple )  )) #find a new NW corner
        push!(corners,(i,j))

        local la=Vector{Int}([])
        local mu=Vector{Int}([])
        local rr=Vector{Int}([])

        local s=0
        while isa( dc.m[i+s,j], Tuple )

          local k=0
          while isa( dc.m[i+s,j+k], Tuple )  # find SE boxes
            k +=1
          end          
          push!(la,k)


          local kk=0
          while j-kk-1>0 && isa( dc.m[i+s,j-kk-1], Tuple )  # find skew boxes
            kk +=1
          end

          mu=mu+fill(kk,length(mu))
          la=la+fill(kk,length(la))
          push!(mu,0)

          if s>0 && i+s>1 && j-kk>1 && ( dc.m[i+s-1,j-kk-1]==1 || dc.m[i+s-1,j-kk-1]==9 || dc.m[i+s-1,j-kk-1]==6  )
            push!(corners,(i+s,j-kk) ) # record new corner
          end
          j=j-kk
          s +=1
        end

        rr=Vector{Vector{Int}}([])
        for el=1:length(la)
          push!(rr, fill(0,mu[el]) )
          for mm=mu[el]+1:la[el]
            push!(rr[end], dc.m[i-1+el,j-1+mm][1]+el )
          end
        end

        push!(lyds,[la,rr,mu,[i-1,j-1]])

      end
    end
  end

  return lyds
end




################
# IN DEVELOPMENT
################


# difference operator (not used)
function ddx(p,i)
# p is a polynomial in x=[x1,x2,..], assumed at least i+1 variables.
  local p1 = evaluate( p, [x[i],x[i+1]], [x[i+1],x[i]] )

  local q=divrem( p-p1, x[i]-x[i+1] )[1]

  return q

end
