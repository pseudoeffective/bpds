# Constructing semistandard tableaux and Schur polynomials
# David Anderson, February 2024.


using Nemo

struct Tableau
  t::Vector{Vector{Int}}
end


# make the row-min superstandard tableau from partition shape la(/mu)
function tableau( la::Vector{Int}, mu=[] )
  tab=[]
  for i in 1:length(la)
    push!(tab, fill(i, la[i]) )
  end

  for i=1:length(mu)
    for j=1:mu[i]
      tab[i][j]=0
    end
  end

  return Tableau(tab)
end


function sstab( la::Vector{Int}, mu=[] )
# make the superstandard tableau from partition shape la(/mu)
  tab=[]
  for i in 1:length(la)
    tr=[]
    if i<=length(mu)
      push!(tr,fill(0,mu[i])...)
      j=mu[i]+1
    else
      j=1
    end
    while j<=la[i]
      k=1
      while k<i && tab[i-k][j]>0
        k+=1
      end
      push!(tr,k)
    j+=1
    end
    push!(tab,tr)
  end

  return Tableau(tab)
end


# extract the (skew) shape of a tableau
function shape(tab::Tableau)
  la=[]
  mu=[]

  len=length(tab.t)

  for i in 1:len
    push!(la,length(tab.t[i]))
    j=1
    while tab.t[i][j]==0
      j+=1
    end
    if j>1
      push!(mu,j-1)
    end  
  end

  return la,mu
end


# add method to display Tableau
function Base.show(io::IO, tab::Tableau)
    println(io)
    for i in 1:length(tab.t)
        for j in 1:length(tab.t[i])
            tt=tab.t[i][j]
            if tt==0
              print(io, 'x', " ")
            else
              print(io, tt, " ")
            end
        end
        println(io)
    end
end



################
# Build tableaux as paths in Young lattice, with memoization

using Memoization

###
# add "a" to row k of tab, return nothing if fails
function add_tabi(tab::Tableau, k, a)
    len = length(tab.t)

    if k > len + 1
        return nothing
    end

    if k == 1
        if a >= last(tab.t[1])
            temp = [[tab.t[1]; a], tab.t[2:end]...]
            return Tableau(temp)
        else
            return nothing
        end
    end

    if k > 1 && k < len + 1
        lam = length(tab.t[k])
        if lam < length(tab.t[k-1]) && 
           a >= last(tab.t[k]) && 
           a > tab.t[k-1][lam+1]
            temp = [tab.t[1:k-1]..., [tab.t[k]; a], tab.t[k+1:end]...]
            return Tableau(temp)
        else
            return nothing
        end
    end

    if k == len + 1
        if a > first(tab.t[len])
            temp = [tab.t..., [a]]
            return Tableau(temp)
        else
            return nothing
        end
    end
end




"""
    ssyt(la::Vector{Int}, ff::Union{Int, Vector{Int}, Vector{Vector{Int}}}; mu::Vector{Int}=[], rowmin::Bool=false) -> Vector{Tableau}

Constructs semistandard Young tableaux (SSYT) based on a given shape and flagging conditions, optionally starting from a skew shape defined by `mu` and with an option for a row-minimal starting point.

## Arguments
- `la::Vector{Int}`: A partition represented as a vector of integers, specifying the desired shape of the tableau. Each integer in the vector represents the number of boxes in each row of the tableau, from top to bottom.
- `ff::Union{Int, Vector{Int}, Vector{Vector{Int}}}`: Flagging conditions for the tableaux. If provided as an integer, it specifies the largest value allowed for an entry. If given as a single vector of integers, it specifies a row flagging condition. If provided as a vector of vectors, it specifies a filling which bounds the tableaux entrywise.
- `mu::Vector{Int}`: An optional vector of integers specifying a subpartition to skew the shape of `la`, representing the rows to be removed from the top of `la` to form a skew tableau. Defaults to an empty vector, indicating no skew.
- `rowmin::Bool`: An optional boolean flag that, when set to `true`, indicates the entries in row `i` should be at least `i`.  Condition is redundant for straight-shape tableaux. Defaults to `false`.

## Returns
- `Vector{Tableau}`: A vector of `Tableau` objects, each representing a semistandard Young tableau that satisfies the given shape, flagging conditions, and optional skew.

# Examples
```julia
### Generate SSYTs for the partition [3, 2] with largest entry 3
ssyts = ssyt([3, 2], 3)

### Generate SSYTs for the shape [3, 2, 1]/[1] and row-minimal condition
ssyts_skewed = ssyt([3, 2, 1], [3, 3, 3], mu=[1], rowmin=true)
"""
@memoize function ssyt(la, bd::Vector{Vector{Int}}; mu = fill(0, length(la)), rowmin=false)
    len = length(la)
    mmu = copy(mu)

    ns = collect( 1:maximum(vcat(bd...)) )
    
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
            push!(tabs, Tableau( [fill(0, mmu[1])] ) )
        end

        if la[1] == mmu[1] + 1
            for i in 1:length(ns)
                if ns[i] <= bd[1][la[1]]
                    push!(tabs, Tableau( [append!( copy(fill(0, mmu[1])), ns[i])] ) )
                end
            end
            return tabs
        end

        if la[1] > mmu[1] + 1
            tabs1 = ssyt([la[1] - 1], bd, mu=mmu)
            for tt in tabs1
                for k in 1:length(ns)
                    if ns[k] <= bd[1][la[1]]
                        t2=add_tabi(tt,1,ns[k])
                        if t2!=nothing
                          push!(tabs, t2 )
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
        
        tabs1 = ssyt(la1, bd, mu=mu1, rowmin=rowmin)
        
        if mmu[end] > 0
            tabs2 = []
            for tt in tabs1
                push!(tabs2, Tableau( vcat( tt.t, [fill(0, mmu[len])] ) ) )
            end
            tabs1 = tabs2
        end
        
        for i in 1:(la[len] - mmu[len])
            tabs2 = []
            for tt in tabs1
                for k in 1:length(ns)
                    if ns[k] <= bd[len][mmu[len]+i] && (!rowmin || ns[k]>=len )
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

###
function ssyt( la, ff::Vector{Int}; mu = fill(0, length(la)), rowmin=false)
# case of flagged tableaux
  bd = Vector{Vector{Int}}([])
  for i=1:length(la)
    push!( bd, fill( ff[i], la[i] ) )
  end

  ssyt( la, bd; mu=mu, rowmin=rowmin )
end


###
function ssyt( la, ff::Int=length(la); mu = fill(0, length(la)), rowmin=false)
# uniformly bounded, default to number of rows of lambda
  bd = Vector{Vector{Int}}([])
  for i=1:length(la)
    push!( bd, fill( ff, la[i] ) )
  end

  ssyt( la, bd; mu=mu, rowmin=rowmin )
end


#####################
# Build tableaux on given shape via iterator, raising operators
#
# Old method via memoization works faster, this is turned off now
#=


function can_increment( tab, bd, i, j )
# determine if i,j entry of tab can be incremented

  # already max
  if tab.t[i][j] >= bd[i][j]
    return false
  end

  # not in shape
  if tab.t[i][j]==0
    return false
  end

  # out of row order
  if j<length(tab.t[i]) && tab.t[i][j] >= tab.t[i][j+1]
    return false
  end

  # out of column order
  if i<length(tab.t) && j <= length(tab.t[i+1]) && tab.t[i][j] >= tab.t[i+1][j]-1
    return false
  end

  return true
end


function increments( tab, bd )
# all tableaux obtained from tab by incrementing a single entry, within bounds bd
  incs = []
  for i in 1:length(tab.t)
    for j in 1:length(tab.t[i])
      if can_increment(tab,bd,i,j)
        tt = deepcopy( tab.t )
        tt[i][j] +=1
        push!(incs, Tableau(tt))
      end
    end
  end
  return incs
end


struct AllSSYT_Iter
  stack::Vector{Any}
  seen::Set{Vector{Vector}}
  bd::Vector{Any}
end

function AllSSYT_Iter( tab::Tableau, bd )
  # initialize
  seen = Set( [tab.t] )
  incs = increments( tab, bd )
  stack = [(tab, incs)]
  return AllSSYT_Iter(stack,seen,bd)
end

Base.IteratorSize(::Type{<:AllSSYT_Iter}) = Base.SizeUnknown()


function Base.iterate(iter::AllSSYT_Iter, state=nothing)

    while !isempty(iter.stack)
        current, incs = pop!(iter.stack)

        unseen_incs = filter( tab -> !(tab.t in iter.seen), incs )

        for tab in unseen_incs
          push!(iter.seen, tab.t)  # mark new tableau as seen
          push!( iter.stack, (tab, increments(tab,iter.bd)) )
        end

        return( current, isempty(iter.stack) ? nothing : iter.stack[end] )
    end

    return nothing  # End of iteration
end


function ssyt_below(tab::Tableau, bd::Vector{Vector{Int}} )
# iterator of all SSYT starting from tab, bounded

  iter = AllSSYT_Iter(tab,bd)

  return iter

end


function ssyt_below(tab::Tableau, ff::Vector{Int} )
# iterator of all SSYT starting from tab, flagged

  bd=Vector{Vector{Int}}([])
  for i in 1:length(tab.t)
    push!(bd, fill(ff[i],length(tab.t[i])))
  end

  return ssyt_below(tab, bd)

end


function ssyt( la::Vector{Int}, ff::Union{Vector{Int},Vector{Vector{Int}}}; mu=[], rowmin=false )
  if rowmin
    tab=tableau( la, mu )
  else
    tab=sstab( la, mu )
  end

  return ssyt_below(tab, ff)
end

=#
################


#################
# Polynomial constructors

######
struct DoublePolyRing
           ring::ZZMPolyRing
           x_vars::Vector{ZZMPolyRingElem}
           y_vars::Vector{ZZMPolyRingElem}
       end

#####


"""
    xy_ring(n,m)

Form a DoublePolyRing object

## Arguments
- `n::Int`: The number of x variables
- `m::Int`: The number of y variables

## Returns
- `DoublePolyRing`, `Vector{ZZMPolyRingElem}`, `Vector{ZZMPolyRingElem}`: The double polynomial ring, with variable sets of specified size.

## Examples

### Define a DoublePolyRing with three x-variables and two y-variables
R,x,y = xy_ring(3,2)

### Extract the x variables

x = R.x_vars

### Extract the y variables

y = R.y_vars

# Notes

The ring is constructed from a `ZZMPolyRing` object in Nemo, and elements can be manipulated as such.

"""
function xy_ring(n,m)
  local xvars = ["x$(i)" for i=1:n]
  local yvars = ["y$(i)" for i=1:m]

  R,all_vars = polynomial_ring(ZZ,vcat(xvars,yvars))

  x = all_vars[1:n]
  y = all_vars[n+1:n+m]

  return DoublePolyRing(R,x,y), x, y

end




# return the product of binomials
function tab2bin( tab::Tableau, RR::DoublePolyRing; xoffset = 0, yoffset = 0 )
  len = length(tab.t)

  x = RR.x_vars
  y = RR.y_vars

  n = length(x)
  m = length(y)

  bin = RR.ring(1)

  for i=1:len
    for j=1:length( tab.t[i] )
      tt = tab.t[i][j]
      if tt>0
        p=RR.ring(0)
        if tt+xoffset<=n
          p=p+x[tt+xoffset]
        end
        if tt+j-i+yoffset<=m && tt+j-i+yoffset>0
          p=p+y[tt+j-i+yoffset]
        end
        bin = bin*p
      end
    end
  end

  return bin

end


# sum of binomials for a set of tableaux
function ssyt2pol( tabs, RR::DoublePolyRing; xoffset=0, yoffset=0 )

  pol=RR.ring(0)

  for tab in tabs
    pol = pol + tab2bin( tab, RR; xoffset=xoffset, yoffset=yoffset )
  end

  return pol

end



"""
    schur_poly(la, ff, RR=xy_ring(length(la), length(la)+la[1])[1]; mu=[], xoffset=0, yoffset=0, rowmin=false)

Compute the Schur polynomial corresponding to a given (skew) partition `la/mu` and a flag `ff`, in an optionally specified ring `RR`. The polynomial is constructed as an enumerator of semistandard Young tableaux of (skew) shape `la/mu` and bounded by the flagging condition `ff`.

## Arguments
- `la::Vector{Int}`: A partition represented as a vector of integers, specifying the shape of the Young diagram.
- `ff::Union{Int,Vector{Int},Vector{Vector{Int}}}`: A flag specifying bounds on the tableaux. If `ff` is given as a single integer, it bounds the entries of the tableaux.  If `ff` is a vector of integers, it must be of length at least that of `la`; then `ff[i]` bounds the entries in the `i`th row of the tableaux.  If `ff` is a vector of vectors, it is interpreted as a tableaux whose shape is assumed to contain `la`; then the entries of `ff` bound the tableaux entrywise.
- `RR::DoublePolyRing`: An optional argument specifying the double polynomial ring to use for constructing the Schur polynomial. Defaults to a ring constructed based on the size of `la`.
- `mu::Vector{Int}`: An optional argument specifying a subpartition of `la`, for skew Schur polynomials. Defaults to an empty vector, for the straight shape `la`.
- `xoffset::Int`: An optional argument specifying an offset value for the x-variable indices in the polynomial. Defaults to 0.
- `yoffset::Int`: An optional argument specifying an offset value for the y-variable indices in the polynomial. Defaults to 0.
- `rowmin::Bool`: An optional argument specifying whether to use row-minimal tableau, i.e., to require that entries in row `i` be at least `i`.  (This is a nontrivial condition only for skew shapes.) Defaults to `false`.

## Returns
- `ZZMPolyRingElem`: The Schur polynomial as an element of the specified polynomial ring `RR`.

# Examples

### Specify a partition
la = [2, 1]

### Specify a flag
ff = 3

### Compute the Schur polynomial
poly = schur_poly(la, ff)

### Display the resulting polynomial
println(poly)

### To get the single Schur polynomial, change the coefficient ring
R,x,y = xy_ring(3,0)

poly = schur_poly(la,ff,R)



# Notes

The function constructs all semistandard Young tableaux that of given shape bounded by the flagging, then sums the corresponding monomials to form the Schur polynomial.

"""
function schur_poly( la, ff::Union{Vector{Int},Vector{Vector{Int}}}, RR::DoublePolyRing=xy_ring( length(la) , length(la)+la[1] )[1]; mu = [], xoffset=0, yoffset=0, rowmin=false )
  if length(la)==0
    return RR.ring(1)
  end

  tbs = ssyt( la, ff, mu=mu, rowmin=rowmin )

  pol = ssyt2pol( tbs, RR; xoffset=xoffset, yoffset=yoffset )

#=
  pol = RR.ring(0)

  for tab in tbs
    pol = pol + tab2bin( tab, RR; xoffset=xoffset, yoffset=yoffset )
  end
=#

  return pol

end



function schur_poly( la, ff::Int, RR::DoublePolyRing=xy_ring( length(la) , length(la)+la[1] )[1]; mu = [], xoffset=0, yoffset=0, rowmin=false )
  if length(la)==0
    return RR.ring(1)
  end

  schur_poly( la, Vector{Int}(fill(ff,length(la))), RR; mu = mu, xoffset=xoffset, yoffset=yoffset, rowmin=rowmin )


end


