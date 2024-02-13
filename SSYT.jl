# Constructing semistandard tableaux and Schur polynomials
# David Anderson, February 2024.


using Nemo

struct Tableau
  t::Vector{Vector{Int}}
end

function tableau( la::Vector{Int}, mu=[] )
# make a tableau from partition shape la(/mu)
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
# make a superstandard tableau from partition shape la(/mu)
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


###
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




#################
# Polynomial constructors

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




function tab2bin( tab::Tableau, RR::DoublePolyRing; xoffset = 0, yoffset = 0 )
# return the product of binomials
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


function ssyt2pol( tabs, RR::DoublePolyRing; xoffset=0, yoffset=0 )
# sum of binomials for a set of tableaux

  pol=RR.ring(0)

  for tab in tabs
    pol = pol + tab2bin( tab, RR; xoffset=xoffset, yoffset=yoffset )
  end

  return pol

end





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


