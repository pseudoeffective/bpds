# Tools for generating BPDs in Julia
# David Anderson, 11 February 2024.



#############
# BPD type and constructors

struct BPD
    m::Matrix{Int8}
end

# Symbol to integer mapping
const SIXVTX_TO_INT = Dict(
    "O" => 0,
    "+" => 1,
    "/" => 2,
    "%" => 3,
    "|" => 4,
    "-" => 5
)

function BPD(matrix::Matrix{String})
    int_matrix = map(x -> SIXVTX_TO_INT[x], matrix)
    return BPD(int_matrix)
end


# convert integers back to symbols for display
function int_to_symbol(i::Int8)
    symbols = ['O', '+', '/', '%', '|', '-']
    return symbols[i+1]  # assuming integers 0-5 map to symbols
end



# add method to Base.show for BPD display
function Base.show(io::IO, bpd::BPD)
    println(io)
    for i in 1:size(bpd.m, 1)
        for j in 1:size(bpd.m, 2)
            print(io, int_to_symbol(bpd.m[i, j]), " ")
        end
        println(io)
    end
end


function Rothe(w)
 # construct the Rothe BPD for w
        local n=length(w)
        local r=Matrix{String}(undef,n,n)

             for j = 1:n
                if j<w[1]
                   r[1,j]="O"
                elseif j==w[1]
                   r[1,j]="/"
                else
                   r[1,j]="-"
                end
             end

             for i=2:n
               for j=1:n
                 if j<w[i]
                   if r[i-1,j]=="/" || r[i-1,j]=="|" || r[i-1,j]=="+"
                      r[i,j]="|"
                   else
                      r[i,j]="O"
                   end
                 elseif j==w[i]
                   r[i,j]="/"
                 else
                   if r[i-1,j]=="|" || r[i-1,j]=="/" || r[i-1,j]=="+"
                     r[i,j]="+"
                   else
                     r[i,j]="-"
                   end
                 end
               end
             end

        return(BPD(r))
end



############
# droop/drip/drop moves

function can_droop(bpd,i1,j1,i2,j2)

 # check rectangle bounds
    if i2<i1+1 || j2<j1+1
       return(false)
    end

 # check NW and SE corners
    if bpd.m[i1,j1] != 2 || bpd.m[i2,j2] != 0
      return(false)
    end

 # check N and S borders
    for j=j1+1:j2
       if bpd.m[i1,j]==2 || bpd.m[i2,j]==2 || bpd.m[i1,j]==3 || bpd.m[i2,j]==3
         return(false)
       end
    end

 # check W and E borders
    for i=i1+1:i2
       if bpd.m[i,j1]==2 || bpd.m[i,j2]==2 || bpd.m[i,j1]==3 || bpd.m[i,j2]==3
         return(false)
       end
    end

 # check inside of rectangle
    for i=i1+1:i2-1
      for j=j1+1:j2-1
        if bpd.m[i,j] == 2 || bpd.m[i,j]==3
          return(false)
        end
      end
    end
  return(true)
end


# to generate bpds by drooping from flats, need to allow larger droops.
# can have 3 on SW or NE corners of rectangle.

function can_flat_drop(bpd,i1,j1,i2,j2)

 # check bounds
    if i2<i1+1 || j2<j1+1
       return(false)
    end

 # small droop not allowed
    if (i2,j2)==(i1+1,j1+1)
       return(false)
    end

 # check corners
    if bpd.m[i1,j1] != 2 || bpd.m[i2,j2] != 0
      return(false)
    end

 # check N and S borders
    for j=j1+1:j2-1
       if bpd.m[i1,j]==2 || bpd.m[i2,j]==2 || bpd.m[i1,j]==3 || bpd.m[i2,j]==3 || bpd.m[i1,j]==0 || bpd.m[i2,j]==0
         return(false)
       end
    end

 # check W and E borders
    for i=i1+1:i2-1
       if bpd.m[i,j1]==2 || bpd.m[i,j2]==2 || bpd.m[i,j1]==3 || bpd.m[i,j2]==3 || bpd.m[i,j1]==0 || bpd.m[i,j2]==0
         return(false)
       end
    end

 # check interior
    for i=i1+1:i2-1
      for j=j1+1:j2-1
        if bpd.m[i,j] == 2 || bpd.m[i,j]==3 || bpd.m[i,j]==0
          return(false)
        end
      end
    end
  return(true)

end


# drip is the small droop
function can_drip(bpd,i1,j1)

 # check corners
    if bpd.m[i1+1,j1+1] != 0
      return(false)
    end

    if bpd.m[i1,j1] == 0 || bpd.m[i1,j1] == 1
      return(false)
    end

    if bpd.m[i1,j1+1] == 0 || bpd.m[i1,j1+1] == 1
      return(false)
    end

    if bpd.m[i1+1,j1] == 0 || bpd.m[i1+1,j1] == 1
      return(false)
    end


  return(true)

end


function droop(bpd,i1,j1,i2,j2)
 # assumes can_[flat_]droop==true
    
    local bpd2=deepcopy(bpd.m)
    # set corners of rectangle
    bpd2[i1,j1]=0
    bpd2[i2,j2]=3
    if bpd2[i1,j2]==5
      bpd2[i1,j2]=2
    elseif bpd2[i1,j2]==3
      bpd2[i1,j2]=4
    end
    if bpd2[i2,j1]==4
      bpd2[i2,j1]=2
    elseif bpd2[i2,j1]==3
      bpd2[i2,j1]=5
    end

    # set west edge
    for i=i1+1:i2-1
       if bpd2[i,j1]==4
          bpd2[i,j1]=0
       elseif bpd2[i,j1]==1
          bpd2[i,j1]=5
       end
    end

    # set north edge
    for j=j1+1:j2-1
       if bpd2[i1,j]==5
          bpd2[i1,j]=0
       elseif bpd2[i1,j]==1
          bpd2[i1,j]=4
       end
    end


    # set east edge
    for i=i1+1:i2-1
       if bpd2[i,j2]==0
          bpd2[i,j2]=4
       elseif bpd2[i,j2]==5
          bpd2[i,j2]=1
       end
    end

    # set south edge
    for j=j1+1:j2-1
       if bpd2[i2,j]==0
          bpd2[i2,j]=5
       elseif bpd2[i2,j]==4
          bpd2[i2,j]=1
       end
    end

    return(BPD(bpd2))
end


function all_droops(bpd)
# produce all droops of bpd
   local n=size(bpd.m)[1]

   local dps = []

   for i1=1:n-1
     for j1=1:n-1
       for i2=i1+1:n
          for j2=j1+1:n
            if can_droop(bpd,i1,j1,i2,j2)
              local bpd2=droop(bpd,i1,j1,i2,j2)
              push!(dps,bpd2)
            end
          end
       end
     end
   end

   return(dps)
end


function flat_drops(bpd)
# produce all (flat) drops of bpd
   local n=size(bpd.m)[1]

   local dps = []

   for i1=1:n-1
     for j1=1:n-1
       for i2=i1+1:n
          for j2=j1+1:n
            if can_flat_drop(bpd,i1,j1,i2,j2)
              bpd2=makeflat(droop(bpd,i1,j1,i2,j2))
              push!(dps,bpd2)
            end
          end
       end
     end
   end

   return(dps)
end


function all_drips(bpd)
# produce all drips of bpd
   local n=size(bpd.m)[1]

   local dps = []

   for i1=1:n-1
     for j1=1:n-1
            if can_drip(bpd,i1,j1)
              local bpd2=droop(bpd,i1,j1,i1+1,j1+1)
              push!(dps,bpd2)
       end
     end
   end

   return(dps)
end


###############
# iterator generating all BPDs for w

struct AllBelowIterator
    stack::Vector{Any}
    seen::Set{Matrix}
end


function AllBelowIterator(bpd::BPD)
    # Initialize with the first element
    seen = Set([bpd.m])
    droops = all_droops(bpd)
    stack = [(bpd, droops)]
    return AllBelowIterator(stack,seen)
end


Base.IteratorSize(::Type{<:AllBelowIterator}) = Base.SizeUnknown()


function Base.iterate(iter::AllBelowIterator, state=nothing)

    while !isempty(iter.stack)
        current, droops = pop!(iter.stack)

        unseen_droops = filter( b -> !(b.m in iter.seen), droops )

        for b in unseen_droops
          push!(iter.seen, b.m)  # mark new droop as seen
          push!( iter.stack, (b, all_droops(b)) )
        end

        return( current, isempty(iter.stack) ? nothing : iter.stack[end] )
    end

    return nothing  # End of iteration
end


function all_bpds(w)
    local bpd = Rothe(w)
    iter = AllBelowIterator(bpd)

    return iter
end


##########


function isflat(bpd)
# determine if bpd is flat
   local n=size(bpd.m)[1]

   for i=2:n-1
     for j=2:n-1
       if bpd.m[i,j]==0 && bpd.m[i-1,j]!=0 && bpd.m[i,j-1]!=0 && bpd.m[i-1,j-1]==2
          return(false)
       end
     end
   end
   return(true)
end


function makeflat(bpd)
# returns flat bpd in drift class of bpd
   local n=size(bpd.m)[1]

   for i=2:n-1
     for j=2:n-1
       if bpd.m[i,j]==0 && bpd.m[i-1,j]!=0 && bpd.m[i,j-1]!=0 && bpd.m[i-1,j-1]==2
         local bpd2=droop(bpd,i-1,j-1,i,j)
         return makeflat(bpd2)
       end
     end
   end

   return(bpd)   

end   


#############
# iterator generating all flat BPDs for w

struct FlatBelowIterator
    stack::Vector{Any}
    seen::Set{Matrix}
end


function FlatBelowIterator(bpd::BPD)
    # Initialize with the first element
    seen = Set([makeflat(bpd).m])
    drops = flat_drops(makeflat(bpd))
    stack = [(makeflat(bpd), drops)]
    return FlatBelowIterator(stack,seen)
end


Base.IteratorSize(::Type{<:FlatBelowIterator}) = Base.SizeUnknown()


function Base.iterate(iter::FlatBelowIterator, state=nothing)

    while !isempty(iter.stack)
        current, drops = pop!(iter.stack)

        unseen_drops = filter( b -> !(makeflat(b).m in iter.seen), drops )

        for b in unseen_drops
          b=makeflat(b)
          push!(iter.seen, b.m)  # mark new drop as seen
          push!( iter.stack, (b, flat_drops(b)) )
        end

        return( makeflat(current), isempty(iter.stack) ? nothing : iter.stack[end] )
    end

    return nothing  # End of iteration
end



function flat_bpds(w)
    local bpd = Rothe(w)
    iter = FlatBelowIterator(bpd)

    return iter
end


function vec_flat_bpds(w)
    local bpd = Rothe(w)
    iter = FlatBelowIterator(bpd)

    return collect(iter)
end



#############
# Conversions to and from ASMs
# where to put these functions?

function bpd2asm( b::BPD )

  local n=size(b.m)[1]

  local a=zeros(Int8,n,n)

  for i=1:n
    for j=1:n
      if b.m[i,j]==2
        a[i,j]=1
      elseif b.m[i,j]==3
        a[i,j]=-1
      end
    end
  end

  return a
end


function asm2bpd( a )

  local n=size(a)[1]

  local b=Matrix{Int8}(undef,n,n)

  if a[n,1]==1
    b[n,1]=2
  else
    b[n,1]=4
  end

  for j=2:n
    if a[n,j]==1
      b[n,j]=2
    elseif b[n,j-1]==2 || b[n,j-1]==1
      b[n,j]=1
    else
      b[n,j]=4
    end
  end

  for i=n-1:-1:1
    for j=1:n
      if a[i,j]==1
        b[i,j]=2

      elseif a[i,j]==-1
        b[i,j]=3

      elseif a[i,j]==0

        if j==1
          if b[i+1,j]==4
            b[i,j]=4
          else
            b[i,j]=0
          end

        else
          if b[i,j-1]==2 || b[i,j-1]==1 || b[i,j-1]==5
            local bb=5
          else
            bb=0
          end
          if b[i+1,j]==3 || b[i+1,j]==1 || b[i+1,j]==4
            local cc=4
          else
            cc=0
          end
          if bb==5
            if cc==4
              b[i,j]=1
            else
              b[i,j]=5
            end
          else
            b[i,j]=cc
          end
        end
      end
    end
  end

  return BPD(b)
end



#########
# old version of iterator, not used
#=

struct AllBelowIterator
    stack::Vector{Any}
    seen::Set{Matrix}
    bpd::BPD
end


function AllBelowIterator(bpd::BPD)
    # Initialize with the first element and an empty set for seen elements
    return AllBelowIterator([(bpd, all_droops(bpd))], Set([bpd.m]), bpd)
end


Base.IteratorSize(::Type{<:AllBelowIterator}) = Base.SizeUnknown()


function Base.iterate(iter::AllBelowIterator, state=nothing)

    if state==nothing
      return( iter.bpd, true )
    end

    while !isempty(iter.stack)
        current, droops = iter.stack[end]

        if isempty(droops)
            pop!(iter.stack)  # Remove the exhausted element
            continue
        end

        next_bpd = popfirst!(droops)  # Get the next droop to process

        if !in(next_bpd.m, iter.seen)
            push!(iter.seen, next_bpd.m)  # Mark as seen
            # Add to stack with its droops
            push!(iter.stack, (next_bpd, all_droops(next_bpd)))
            return next_bpd, iter.stack
        end
    end

    return nothing  # End of iteration
end


function all_bpds(w)
    local bpd = Rothe(w)
    iter = AllBelowIterator(bpd)

    return iter
end

=#
#########

####################
#not used
function can_sharp_drop(bpd,i1,j1,i2,j2)

 # check bounds
    if i2<i1+1 || j2<j1+1
       return(false)
    end

 # small droop not allowed
    if (i2,j2)==(i1+1,j1+1)
       return(false)
    end

 # check corners
    if bpd.m[i1,j1] != 2 || bpd.m[i2,j2] != 0
      return(false)
    end

 # check active pipe
    if i2==i1+1 && bpd.m[i2-1,j2-1] == 5
      return(false)
    end
    if j2==j1+1 && bpd.m[i2-1,j2-1]==4
      return(false)
    end

 # check NW of destination
    if bpd.m[i2-1,j2-1]==0
      return(false)
    end

 # check N and S borders
    for j=j1+1:j2-1
       if bpd.m[i1,j]==2 || bpd.m[i2,j]==2 || bpd.m[i1,j]==3 || bpd.m[i2,j]==3
         return(false)
       end
    end

 # check W and E borders
    for i=i1+1:i2-1
       if bpd.m[i,j1]==2 || bpd.m[i,j2]==2 || bpd.m[i,j1]==3 || bpd.m[i,j2]==3
         return(false)
       end
    end

 # check interior
    for i=i1+1:i2-1
      for j=j1+1:j2-1
        if bpd.m[i,j] == 2 || bpd.m[i,j]==3
          return(false)
        end
      end
    end
  return(true)
end



function sharp_drops(bpd)
# produce all (sharp) drops of bpd
   local n=size(bpd.m)[1]

   local dps = []

   for i1=1:n-1
     for j1=1:n-1
       for i2=i1+1:n
          for j2=j1+1:n
            if can_sharp_drop(bpd,i1,j1,i2,j2)
              bpd2=droop(bpd,i1,j1,i2,j2)
              push!(dps,bpd2)
            end
          end
       end
     end
   end

   return(dps)
end



function issharp(bpd)
# determine if bpd is sharp
   local n=size(bpd.m)[1]

   for i=1:n-1
     for j=1:n-1
       if bpd.m[i,j]==0 && bpd.m[i+1,j]!=0 && bpd.m[i,j+1]!=0 && bpd.m[i+1,j+1]!=1
          return(false)
       end
     end
   end
   return(true)
end



function makesharp(bpd)
# returns sharp bpd in drift class of bpd
   local n=size(bpd.m)[1]

   for i=1:n-2
     for j=1:n-2
       if bpd.m[i,j]==0 && bpd.m[i+1,j]!=0 && bpd.m[i,j+1]!=0 && bpd.m[i+1,j+1]==3
         local bpd2=deepcopy(bpd.m)
         bpd2[i,j]=2
         bpd2[i+1,j+1]=0

         if bpd2[i+1,j]==2
           bpd2[i+1,j]=4
         else
           bpd2[i+1,j]=3
         end

         if bpd2[i,j+1]==2
           bpd2[i,j+1]=5
         else
           bpd2[i,j+1]=3
         end

         return makesharp(BPD(bpd2))
       end
     end
   end

   return(bpd)   

end

