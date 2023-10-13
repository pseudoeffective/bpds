# Tools for computing Schubert polynomials in Julia, using BPDs
# David Anderson, October 2023.


#Using Nemo

include("ssyt.jl")
include("BpdBase.jl")

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


function bpd2bin( bpd, R::DoublePolyRing=xy_ring( size(bpd)[1]-1, size(bpd)[2]-1 )[1]  )
# product of binomials for bpd
# requires DoublePolyRing
# can get single polyn by using no y_vars
  local n=size(bpd)[1]-1
  bin = R.ring(1)

  x = R.x_vars
  y = R.y_vars

  local aa=length(x)
  local bb=length(y)

  for i=1:n
    for j=1:n
      if bpd[i,j]=="O"
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



function bpd2sd( bpd, R::DoublePolyRing=xy_ring( size(bpd)[1]-1, size(bpd)[2]-1 )[1]  )
# compute drift class polynomial from flat bpd
  if !isflat(bpd)
    return 0
  end

  sd = R.ring(1)

  tcomps = tableau_components(bpd)

  for tt in tcomps
    tt[2] = tt[2]+collect(1:length(tt[2]))
    sd = sd*schur_poly( tt[1], tt[2], R; mu=tt[3], xoffset=tt[4][1], yoffset=tt[4][2], rowmin=true )
  end

  return sd

end


function schub_flats( w, R::DoublePolyRing=xy_ring( length(w)-1, length(w)-1 )[1] )
# drift class formula

  fbpds = flat_bpds(w)

  pol = R.ring(0)

  for b in fbpds
    pol = pol+bpd2sd(b,R)
  end

  return pol

end

