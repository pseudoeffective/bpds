# Tools for computing Schubert polynomials in Julia, using BPDs
# David Anderson, October 2023.


#Using Nemo

#include("ssyt.jl")
include("SSYT.jl")
include("BpdBase.jl")
include("Drifts.jl")



"""
    schub_poly(w, R; method)

Return the Schubert polynomial for the permutation `w`

## Arguments
- `w::Vector{Int}`: A permutation
- `R::DoublePolyRing`: The ambient polynomial ring, with underlying ring `R.ring` and variables `R.x_vars` and `R.y_vars`
- `method`: An optional argument specifying the algorithm for computing.

The options for `method` are:
- `method="drift"` (default), computes from flat BPDs by drift class formula
- `method="bpd"`, computes by summing over all BPDs
- `method="dd"`, computes by divided differences from longest permutation 

## Returns
`ZZMPolyRingElem`: the Schubert polynomial in the ring R

## Examples
```julia
# Define a single-variable DoublePolyRing
R,x,y = xy_ring(5,0)

# Choose a permutation
w = [3,5,1,4,2]

# Compute via drift classes (default)
pol1 = schub_poly(w, R)

# Compute via divided differences
pol2 = schub_poly(w, R, method="dd")

pol1==pol2
```

### If there are not enough x-variables to compute by descending induction, the "dd" method throws an error
```julia
R,x,y = xy_ring(3,0);

schub_poly(w, R, method="dd");

```

### If R is not specified, a double polynomial ring containing the double Schubert polynomial is chosen

```julia
pol3 = schub_poly(w, method="bpd");

R,x,y = xy_ring(4,4);

pol4 = schub_poly(w, R, method="drift");

pol3==pol4
```
"""
function schub_poly(w, R::DoublePolyRing=xy_ring( max(length(w)-1,1), max(length(w)-1,1) )[1]; method="drift" )

  if method=="dd"
    return schub_dd(w,R)
  elseif method=="bpd"
    return schub_bpd(w,R)
  end

  return schub_drifts(w,R)

end


## TO DO: Add "transition" method for schub_poly


"""
    groth_poly(w, R; method)

Return the Grothendieck polynomial for the permutation `w`

## Arguments
- `w::Vector{Int}`: A permutation
- `R::DoublePolyRing`: The ambient polynomial ring, with underlying ring `R.ring` and variables `R.x_vars` and `R.y_vars`
- `method`: An optional argument specifying the algorithm for computing.

The options for `method` are:
- `method="bpd"` (default), computes by summing over all BPDs
- `method="dd"`, computes by divided differences from longest permutation 

## Returns
`ZZMPolyRingElem`: the Grothendieck polynomial in the ring R

## Examples
```julia
# Define a single-variable DoublePolyRing
R,x,y = xy_ring(5,0)

# Choose a permutation
w = [3,5,1,4,2]

# Compute via bpds (default)
pol1 = groth_poly(w, R)

# Compute via divided differences
pol2 = groth_poly(w, R, method="dd")

pol1==pol2
```
### If R is not specified, a double polynomial ring containing the double Schubert polynomial is chosen

```julia
pol3 = groth_poly(w, method="bpd");

R,x,y = xy_ring(4,4);

pol4 = groth_poly(w, R, method="dd");

pol3==pol4
```
"""
function groth_poly(w, R::DoublePolyRing=xy_ring( max(length(w)-1,1), max(length(w)-1,1) )[1]; method="bpd" )

  if method=="dd"
    return groth_dd(w,R)
  elseif method=="drift"
    return groth_drifts(w,R)
  end

  return groth_bpd(w,R)

end


## TO DO: add "transition" method for groth_poly
## TO DO: figure out "drift" method for groth_poly



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


function bpd2bin( bpd, R::DoublePolyRing=xy_ring( size(bpd.m)[1]-1, size(bpd.m)[2]-1 )[1]; version="schub"  )
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
        if version=="groth" && i<=aa && j<=bb
          p=p-x[i]*y[j]
        end     
        bin = bin*p
        if version=="groth"
          bin = -bin
        end
      end

      if version=="groth" && bpd.m[i,j]==3
        p=R.ring(1)
        if i<=aa
          p=p-x[i]
        end
        if j<=bb
          p=p-y[j]
        end
        if i<=aa && j<=bb
          p=p+x[i]*y[j]
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


function schub_bpd( w, R::DoublePolyRing=xy_ring( max(length(w)-1,1), max(length(w)-1,1) )[1]  )
# compute schubert pol by bpd formula
  bpds=all_bpds(w)

  pol=R.ring(0)

  for bp in bpds
    pol = pol+bpd2bin(bp,R)
  end

  return(pol)

end




function schub_drifts( w, R::DoublePolyRing=xy_ring( max(length(w)-1,1), max(length(w)-1,1) )[1] )
# compute schubert pol by drift class formula

  fbpds = flat_bpds(w)

  pol = R.ring(0)

  for b in fbpds
    b=markconfig(b)
    pol = pol+dc2sd(b,R)
  end

  return pol

end


@memoize function schub_dd( w, R::DoublePolyRing=xy_ring( max(length(w)-1,1), max(length(w)-1,1) )[1] )
# compute schubert pol by divided differences

  n=length(w)
  while n>0 && w[n]==n
    pop!(w)
    n=length(w)
  end

  if n==0
    return(R.ring(1))
  end

  if length(R.x_vars)<n-1
    throw(ArgumentError("Not enough x variables for this method"))
  end

  i=n-1
  while i>0 && w[i]>w[i+1]
    i=i-1
  end

  if i==0
    return sp0( n, R )
  end

  w = vcat( w[1:i-1], w[i+1], w[i], w[i+2:n] )

  pol1 = schub_dd( w, R )

  return ddx( pol1, i, R )

end




function groth_bpd( w, R::DoublePolyRing=xy_ring( max(length(w)-1,1), max(length(w)-1,1) )[1]  )
# compute grothendieck pol by bpd formula
  bpds=all_Kbpds(w)

  pol=R.ring(0)

  for bp in bpds
    pol = pol+bpd2bin(bp,R, version="groth")
  end

  return((-1)^len(w)*pol)

end


@memoize function groth_dd( w, R::DoublePolyRing=xy_ring( max(length(w)-1,1), max(length(w)-1,1) )[1] )
# compute grothendieck pol by divided differences

  n=length(w)
  while n>0 && w[n]==n
    pop!(w)
    n=length(w)
  end

  if n==0
    return(R.ring(1))
  end

  if length(R.x_vars)<n-1
    throw(ArgumentError("Not enough x variables for this method"))
  end

  i=n-1
  while i>0 && w[i]>w[i+1]
    i=i-1
  end

  if i==0
    return gp0( n, R )
  end

  w = vcat( w[1:i-1], w[i+1], w[i], w[i+2:n] )

  pol1 = groth_dd( w, R )

  return pdx( pol1, i, R )

end





#############

function sp0( n, R::DoublePolyRing=xy_ring( n-1, n-1 )[1] )

  x=R.x_vars
  y=R.y_vars

  aa=length(x)
  bb=length(y)

  pol=1

  for i in 1:n-1
    for j in 1:n-i
      p=R.ring(0)
      if i<=aa p=p+x[i] end
      if j<=bb p=p+y[j] end
      pol = pol*p
    end
  end

  return(pol)

end


function gp0( n, R::DoublePolyRing=xy_ring( n-1, n-1 )[1] )

  x=R.x_vars
  y=R.y_vars

  aa=length(x)
  bb=length(y)

  pol=1

  for i in 1:n-1
    for j in 1:n-i
      p=R.ring(0)
      if i<=aa p=p+x[i] end
      if j<=bb p=p+y[j] end
      if i<=aa && j<=bb p=p-x[i]*y[j] end
      pol = pol*p
    end
  end

  return(pol)

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




# difference operator
function ddx(p,i,R)
# p is a polynomial in R.x_vars

  x = R.x_vars

  if i>length(x)
    return(R.ring(0))
  end

  if i==length(x)
    p1 = evaluate( p, [x[i]], [R.ring(0)] )
    q=divrem( p-p1, x[i] )[1]
    return q
  end

  p1 = evaluate( p, [x[i],x[i+1]], [x[i+1],x[i]] )

  q=divrem( p-p1, x[i]-x[i+1] )[1]

  return q

end


# isobaric difference operator
function pdx(p,i,R)
# p is a polynomial in R.x_vars

  x = R.x_vars

  if i>length(x)
    return(p)
  end

  if i==length(x)
    p1 = evaluate( p, [x[i]], [R.ring(0)] )
    q=divrem( p-(1-x[i])*p1, x[i] )[1]
    return q
  end

  p1 = evaluate( p, [x[i],x[i+1]], [x[i+1],x[i]] )

  q=divrem( (1-x[i+1])*p-(1-x[i])*p1, x[i]-x[i+1] )[1]

  return q

end