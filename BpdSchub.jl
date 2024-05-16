# Tools for computing Schubert polynomials in Julia, using BPDs
# David Anderson, April 2024


#Using Nemo
Using LinearAlgebra

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
- `memo`: An optional argument for memoization (default=false)

The options for `method` are:
- `method="transition"`, computes via transition formula
- `method="drift"`, computes from flat BPDs by drift class formula
- `method="bpd"`, computes by summing over all BPDs
- `method="dd"`, computes by divided differences from longest permutation 
- `method="auto" (default), computes by transition if there are y variables, or if w has relatively small length, otherwise computes by divided difference

## Returns
`ZZMPolyRingElem`: the Schubert polynomial in the ring R

## Examples
```julia
# Define a single-variable DoublePolyRing
julia> R = xy_ring(5,0)[1];

# Choose a permutation
julia> w = [3,5,1,4,2];

# Compute via transition classes
julia> pol1 = schub_poly(w, R, method="transition");

# Compute via divided differences
julia> pol2 = schub_poly(w, R, method="dd");

julia> pol1==pol2
true
```

### If there are not enough x-variables to compute by descending induction, the "dd" method throws an error
```julia
julia> R = xy_ring(3,0)[1];

julia> schub_poly(w, R, method="dd")
ERROR: ArgumentError: Not enough x variables for this method

```

### If R is not specified, a double polynomial ring containing the double Schubert polynomial is chosen

```julia
julia> pol3 = schub_poly(w, method="bpd");

julia> R,x,y = xy_ring(4,4);

julia> pol4 = schub_poly(w, R, method="drift");

julia> pol3==pol4
true
```
"""
function schub_poly(w, R::DoublePolyRing=xy_ring( max(length(w)-1,1), max(length(w)-1,1) )[1]; method="auto", memo=false )

  w=trimw(w)
  ws = cutw(w)

  if length(ws)>1
    return( schub_poly(ws[1], R; method=method, memo=memo )*schub_poly(ws[2], R; method=method, memo=memo ) ) 
  end

  if memo
    return schub_trans_memo(w,R)
  end

  if method=="dd"
    return schub_dd(w,R)
  elseif method=="bpd"
    return schub_bpd(w,R)
  elseif method=="drift"
    return schub_drifts(w,R)
  elseif method=="transition"
    return schub_trans(w,R)
  end

  if length(R.y_vars)>0 || len(w)<.2*length(w)^2
    return schub_trans(w,R)
  end

  return schub_dd(w,R)

end


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
julia> R,x,y = xy_ring(5,0);

# Choose a permutation
julia> w = [3,5,1,4,2];

# Compute via bpds (default)
julia> pol1 = groth_poly(w, R);

# Compute via divided differences
julia> pol2 = groth_poly(w, R, method="dd");

julia> pol1==pol2
true
```
### If R is not specified, a double polynomial ring containing the double Schubert polynomial is chosen

```julia
julia> pol3 = groth_poly(w, method="bpd");

julia> R,x,y = xy_ring(4,4);

julia> pol4 = groth_poly(w, R, method="dd");

julia> pol3==pol4
true
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


######
# just counting terms



function nschub(w)
# count the number of terms in the Schubert polynomial
  w = trimw(w)

  if length(w)==0
    return 1
  end

  ws = cutw(w)

  if length(ws)>1
    return( nschub(ws[1])*nschub(ws[2]) ) 
  end

  mtx = max_transition(w)

  if mtx==[]
    b=makeflat(Rothe(w))
    b=markconfig(b)
    tt=tableau_components(b)[1]

    ff = Int[]
    for i in 1:length(tt[2])
      push!(ff, last(tt[2][i]) )
    end

    return round(Int, vex_det( tt[1], ff ) )
  end

  vv=mtx[2]
  wws=mtx[3]

  ss = nschub(vv)

  for ww in wws
    ss = ss+nschub(ww)
  end

  return ss
end



function ngroth(w)
# count the number of terms in the Grothendieck polynomial
  w = trimw(w)

  n = length(w)-1

  R = xy_ring(n,0)[1]

  pw = groth_poly(w,R,method="dd")

  return sum(abs.(coefficients(pw)))  

end

######


function bpd2bin( bpd::BPD, R::DoublePolyRing=xy_ring( size(bpd.m)[1]-1, size(bpd.m)[2]-1 )[1]; version="schub"  )
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


#@memoize 
function schub_trans( w, R::DoublePolyRing=xy_ring( max(length(w)-1,1), max(length(w)-1,1) )[1] )
# compute schubert pol by transition formula

  mxt = max_transition(w)

  if length(mxt)==0
    return schub_drifts(w,R)
  end

  (r,s) = mxt[1]
  vv = mxt[2]
  wws = mxt[3]

  x = R.x_vars
  y = R.y_vars

  pv = schub_trans(vv, R)

  p1 = R.ring(0)

  if r<=length(x)
    p1 = p1+x[r]
  end

  if w[s]<=length(y)
    p1 = p1 + y[w[s]]
  end

  p1 = p1*pv

  for ww in wws
    pw = schub_trans(ww,R)
    p1 = p1+pw
  end

  return p1

end

@memoize function schub_trans_memo( w, R::DoublePolyRing=xy_ring( max(length(w)-1,1), max(length(w)-1,1) )[1] )
# compute schubert pol by transition formula, memoized

  mxt = max_transition(w)

  if length(mxt)==0
    return schub_drifts(w,R)
  end

  (r,s) = mxt[1]
  vv = mxt[2]
  wws = mxt[3]

  x = R.x_vars
  y = R.y_vars

  pv = schub_trans_memo(vv, R)

  p1 = R.ring(0)

  if r<=length(x)
    p1 = p1+x[r]
  end

  if w[s]<=length(y)
    p1 = p1 + y[w[s]]
  end

  p1 = p1*pv

  for ww in wws
    pw = schub_trans_memo(ww,R)
    p1 = p1+pw
  end

  return p1

end


function schub_dd( w, R::DoublePolyRing=xy_ring( max(length(w)-1,1), max(length(w)-1,1) )[1] )
# compute schubert pol by divided differences

  w=trimw(w)

  n=length(w)

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
# drift configuration to s-polynomial
# must take marked configuration as input

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
# return labelled tableaux for a drift configuration dc
# must take marked configuration as input

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



#####################
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



###################
# permutation tools


function len( w::Vector{Int} )
# coxeter length of a permutation (or word)
  n = length(w)
  a = 0

  for i in 1:n-1
    for j in i+1:n
      if w[i]>w[j]
        a=a+1
      end
    end
  end

  return a
end


function cutw( w::Vector{Int} )
# split a permutation w = w1 w2 into two which permute distinct indices, if possible

  n=length(w)

  if n<3
    return [w]
  end

  for i in 1:n-3
    if Set(w[n-i:n])==Set(n-i:n)
      w1 = w[1:n-i-1]
      w2 = vcat(1:n-i-1,w[n-i:n])
      if w1==collect(1:n-i-1)
        return [w2]
      end
      return [w1,w2]
    end
  end

  return [w]

end


function findlast132( w )
# find the lex-last subsequence of pattern 132, return its indices or [] if none

  n=length(w)

  if n<3
    return []
  end

  if n==3
    w==[1,3,2] && return [1,2,3]
    return []
  end

  for i in n-2:-1:1
    for j in n-1:-1:i+1
      for k in n:-1:j+1
        if invperm(sortperm( w[[i,j,k]] ) )==[1,3,2]
          return [i,j,k]
        end
      end
    end
  end

  return []

end


function findlast2143( w )
# find the lex-last subsequence of pattern 2143, return its indices or [] if none

  n=length(w)

  if n<4
    return []
  end

  if n==4
    w==[2,1,4,3] && return [1,2,3,4]
    return []
  end

  for i in n-3:-1:1
    for j in n-2:-1:i+1
      for k in n-1:-1:j+1
        for l in n:-1:k+1
          if invperm(sortperm( w[[i,j,k,l]] ) )==[2,1,4,3]
            return [i,j,k,l]
          end
        end
      end
    end
  end

  return []

end


function maxinversion(w)

  n=length(w)

  if n<2 return [] end
  if n==2
    w==[2,1] && return [1,2]
    return []
  end

  for i in n-1:-1:1
    for j in n:-1:i+1
      if w[i]>w[j]
        return [i,j]
      end
    end
  end

  return []

end



function sij( w, i, j )
# transposition w t_ij

  n=length(w)
  if i==j return w end

  if i>j
    a=j
    b=i
    return sij(w,a,b)
  end

  wt = vcat( w[1:i-1], w[j], w[i+1:j-1], w[i], w[j+1:n] )

  return wt

end



function dominant_transition( w )
# do step in Fan-Guo-Sun transition tree

  pos = findlast132(w)

  if length(pos)==0
    return []
  end

  r,s = pos[2],pos[3]

  vv = sij(w,r,s)

  ws = []
  ll = len(w)

  for i in 1:r-1
    ww = sij(vv, i, r )
    if len(ww)==ll
      push!(ws,ww)
    end
  end

  return [(r,s),vv, ws]

end


function max_transition( w )
# do step in Lascoux-Schutzenberger transition tree

  pos = findlast2143(w)

  if length(pos)==0
    return []
  end

  r,s = pos[3],pos[4]

  vv = sij(w,r,s)

  ws = []
  ll = len(w)

  for i in 1:r-1
    ww = sij(vv, i, r )
    if len(ww)==ll
      push!(ws,ww)
    end
  end

  return [(r,s),vv, ws]

end


function trimw( w )
# remove trailing w(i)=i
  n=length(w)

  ww = deepcopy(w)
  while n>0 && ww[n]==n
    pop!(ww)
    n=length(ww)
  end

  return ww

end

# for counting vexillary terms
function vex_det( la, ff )

  n=length(la)

  if n==0 return 1 end

#  A = zeros(BigFloat,n,n) # may need this depending on float tolerance
  A = zeros(n,n)

  for i in 1:n
    for j in 1:n
      if la[i]-i+j>=0
        A[i,j]=binomial( la[i]-i+j+ff[i]-1, ff[i]-1 )
      end
    end
  end

  return det(A)

end
