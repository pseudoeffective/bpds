# Tools for generating BPDs in Julia
# David Anderson, October 2023.


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

        return(r)
end



function can_droop(bpd,i1,j1,i2,j2)

 # check rectangle bounds
    if i2<i1+1 || j2<j1+1
       return(false)
    end

 # check NW and SE corners
    if bpd[i1,j1] != "/" || bpd[i2,j2] != "O"
      return(false)
    end

 # check N and S borders
    for j=j1+1:j2
       if bpd[i1,j]=="/" || bpd[i2,j]=="/" || bpd[i1,j]=="%" || bpd[i2,j]=="%"
         return(false)
       end
    end

 # check W and E borders
    for i=i1+1:i2
       if bpd[i,j1]=="/" || bpd[i,j2]=="/" || bpd[i,j1]=="%" || bpd[i,j2]=="%"
         return(false)
       end
    end

 # check inside of rectangle
    for i=i1+1:i2-1
      for j=j1+1:j2-1
        if bpd[i,j] == "/" || bpd[i,j]=="%"
          return(false)
        end
      end
    end
  return(true)
end



# to generate bpds by drooping from flats, need to allow larger droops.
# can have "%" on SW or NE corners of rectangle.

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
    if bpd[i1,j1] != "/" || bpd[i2,j2] != "O"
      return(false)
    end

 # check N and S borders
    for j=j1+1:j2-1
       if bpd[i1,j]=="/" || bpd[i2,j]=="/" || bpd[i1,j]=="%" || bpd[i2,j]=="%" || bpd[i1,j]=="O" || bpd[i2,j]=="O"
         return(false)
       end
    end

 # check W and E borders
    for i=i1+1:i2-1
       if bpd[i,j1]=="/" || bpd[i,j2]=="/" || bpd[i,j1]=="%" || bpd[i,j2]=="%" || bpd[i,j1]=="O" || bpd[i,j2]=="O"
         return(false)
       end
    end

 # check interior
    for i=i1+1:i2-1
      for j=j1+1:j2-1
        if bpd[i,j] == "/" || bpd[i,j]=="%" || bpd[i,j]=="O"
          return(false)
        end
      end
    end
  return(true)

end


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
    if bpd[i1,j1] != "/" || bpd[i2,j2] != "O"
      return(false)
    end

 # check active pipe
    if i2==i1+1 && bpd[i2-1,j2-1] == "-"
      return(false)
    end
    if j2==j1+1 && bpd[i2-1,j2-1]=="|"
      return(false)
    end

 # check NW of destination
    if bpd[i2-1,j2-1]=="O"
      return(false)
    end

 # check N and S borders
    for j=j1+1:j2-1
       if bpd[i1,j]=="/" || bpd[i2,j]=="/" || bpd[i1,j]=="%" || bpd[i2,j]=="%"
         return(false)
       end
    end

 # check W and E borders
    for i=i1+1:i2-1
       if bpd[i,j1]=="/" || bpd[i,j2]=="/" || bpd[i,j1]=="%" || bpd[i,j2]=="%"
         return(false)
       end
    end

 # check interior
    for i=i1+1:i2-1
      for j=j1+1:j2-1
        if bpd[i,j] == "/" || bpd[i,j]=="%"
          return(false)
        end
      end
    end
  return(true)
end


function droop(bpd,i1,j1,i2,j2)
 # assumes can_[flat_]droop==true
    
    local bpd2=deepcopy(bpd)
    # set corners of rectangle
    bpd2[i1,j1]="O"
    bpd2[i2,j2]="%"
    if bpd2[i1,j2]=="-"
      bpd2[i1,j2]="/"
    elseif bpd2[i1,j2]=="%"
      bpd2[i1,j2]="|"
    end
    if bpd2[i2,j1]=="|"
      bpd2[i2,j1]="/"
    elseif bpd2[i2,j1]=="%"
      bpd2[i2,j1]="-"
    end

    # set west edge
    for i=i1+1:i2-1
       if bpd2[i,j1]=="|"
          bpd2[i,j1]="O"
       elseif bpd2[i,j1]=="+"
          bpd2[i,j1]="-"
       end
    end

    # set north edge
    for j=j1+1:j2-1
       if bpd2[i1,j]=="-"
          bpd2[i1,j]="O"
       elseif bpd2[i1,j]=="+"
          bpd2[i1,j]="|"
       end
    end


    # set east edge
    for i=i1+1:i2-1
       if bpd2[i,j2]=="O"
          bpd2[i,j2]="|"
       elseif bpd2[i,j2]=="-"
          bpd2[i,j2]="+"
       end
    end

    # set south edge
    for j=j1+1:j2-1
       if bpd2[i2,j]=="O"
          bpd2[i2,j]="-"
       elseif bpd2[i2,j]=="|"
          bpd2[i2,j]="+"
       end
    end

    return(bpd2)
end



function all_droops(bpd)
# produce all droops of bpd
   local n=size(bpd)[1]

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
   local n=size(bpd)[1]

   local dps = Vector{Matrix}([])

   for i1=1:n-1
     for j1=1:n-1
       for i2=i1+1:n
          for j2=j1+1:n
            if can_flat_drop(bpd,i1,j1,i2,j2)
              bpd2=flatten(droop(bpd,i1,j1,i2,j2))
              push!(dps,bpd2)
            end
          end
       end
     end
   end

   return(dps)
end


function sharp_drops(bpd)
# produce all (sharp) drops of bpd
   local n=size(bpd)[1]

   local dps = Vector{Matrix}([])

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




# hash table for all_below
hash_all_below = Dict{Matrix, Set}()

function all_below(bpd)
# returns set of all diagrams below bpd in droop order

   if haskey(hash_all_below, bpd)
     return hash_all_below[bpd]
   end

   local dps=Set( all_droops(bpd) )

   if length(dps)==0
      hash_all_below[bpd]=Set( [bpd] )
      return hash_all_below[bpd]
   end

   local alldps = Set([])

   for b in dps
     alldps=union(alldps,all_below(b))
   end

   push!(alldps,bpd)

   hash_all_below[bpd]=alldps

   return alldps
end



function all_bpds(w)
# returns vector of all bpds for w

  local bpd=Rothe(w)

  bs = all_below(bpd)

  empty!(hash_all_below)  # clear the lookup table

  return collect( bs )

end


function issharp(bpd)
# determine if bpd is sharp
   local n=size(bpd)[1]

   for i=1:n-1
     for j=1:n-1
       if bpd[i,j]=="O" && bpd[i+1,j]!="O" && bpd[i,j+1]!="O" && bpd[i+1,j+1]!="+"
          return(false)
       end
     end
   end
   return(true)
end

function isflat(bpd)
# determine if bpd is flat
   local n=size(bpd)[1]

   for i=2:n-1
     for j=2:n-1
       if bpd[i,j]=="O" && bpd[i-1,j]!="O" && bpd[i,j-1]!="O" && bpd[i-1,j-1]=="/"
          return(false)
       end
     end
   end
   return(true)
end


function flatten(bpd)
# returns flat bpd in drift class of bpd
   local n=size(bpd)[1]

   for i=2:n-1
     for j=2:n-1
       if bpd[i,j]=="O" && bpd[i-1,j]!="O" && bpd[i,j-1]!="O" && bpd[i-1,j-1]=="/"
         local bpd2=droop(bpd,i-1,j-1,i,j)
         return flatten(bpd2)
       end
     end
   end

   return(bpd)   

end   


function sharpen(bpd)
# returns sharp bpd in drift class of bpd
   local n=size(bpd)[1]

   for i=1:n-2
     for j=1:n-2
       if bpd[i,j]=="O" && bpd[i+1,j]!="O" && bpd[i,j+1]!="O" && bpd[i+1,j+1]=="%"
         local bpd2=deepcopy(bpd)
         bpd2[i,j]="/"
         bpd2[i+1,j+1]="O"

         if bpd2[i+1,j]=="/"
           bpd2[i+1,j]="|"
         else
           bpd2[i+1,j]="%"
         end

         if bpd2[i,j+1]=="/"
           bpd2[i,j+1]="-"
         else
           bpd2[i,j+1]="%"
         end

         return sharpen(bpd2)
       end
     end
   end

   return(bpd)   

end



function flat_below(bpd)
# returns vector of all flat bpds below bpd in drop order
# assumes bpd is flat

   local dps=flat_drops(bpd)

   if length(dps)==0
     return [bpd]
   end

   local alldps = Vector{Matrix}([])

   for b in dps
     alldps=union(alldps,flat_below(b))
   end

   prepend!(alldps,[bpd])

   return alldps
end



function flat_bpds(w)
# returns all flat bpds for w

  local bpd=flatten(Rothe(w))

  return flat_below(bpd)

end





################
# IN DEVELOPMENT
################


function drift_class(w,bpd)
# return vector of all diagrams for w in drift class of bpd
# currently very inefficient
# TO DO: improve by writing a drip function, which does the drift moves

  local bflat=flatten(bpd)

  dr = []
  
  for b in all_bpds(w)
    if flatten(b)==bflat
      push!(dr,b)
    end
  end

  return dr
end




# possibly move this one to bpd-schub

function tableau_components(bpd)
# return labelled tableaux for a flat bpd

  local n=size(bpd)[1]

  if !isflat(bpd)
    return( tableau_components( flatten(bpd) ) )
  end

  local lyds=Vector{Vector}([])

  local corners=Vector{Tuple{Int,Int}}([])

  for i=1:n
    for j=1:n
      if !( (i,j) in corners) && bpd[i,j]=="O" && ((i,j)==(1,1) || (i>1 && j>1 &&bpd[i-1,j-1]=="+")) #find a new NW corner
        push!(corners,(i,j))

        local la=Vector{Int}([])
        local mu=Vector{Int}([])
        local rr=Vector{Int}([])

        local s=0
        while bpd[i+s,j]=="O"

          local k=0
          while bpd[i+s,j+k]=="O"  # find SE boxes
            k +=1
          end          
          push!(la,k)

          local el=1
            while bpd[i+s+el,j+k-1+el]=="%" || bpd[i+s+el,j+k-1+el]=="|"
              el +=1
            end
          push!(rr,el-1)


          local kk=0
          while j-kk-1>0 && bpd[i+s,j-kk-1]=="O"  # find skew boxes
            kk +=1
          end

          mu=mu+fill(kk,length(mu))
          la=la+fill(kk,length(la))
          push!(mu,0)

          if s>0 && i+s>1 && j-kk>1 && bpd[i+s-1,j-kk-1]=="+"
            push!(corners,(i+s,j-kk) ) # record new corner
          end
          j=j-kk
          s +=1
        end

#=
        for k=length(la):-1:2
          if la[k-1]==la[k]
            rr[k-1]=rr[k]
          end
        end
=#

        push!(lyds,[la,rr,mu,[i-1,j-1]])

      end
    end
  end

  return lyds
end


#############
# for testing
#############
function sharps(w)

  filter(issharp, all_bpds(w))

end

function flats(w)

  filter(isflat, all_bpds(w))

end

