using Memoization

function Rothe(w)
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

function can_flat_droop(bpd,i1,j1,i2,j2)

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

 # OMIT
 # check destination is a corner
 #   if bpd[i2+1,j2]=="O" || bpd[i2,j2+1]=="O"
 #      return false
 #   end

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


function can_sharp_droop(bpd,i1,j1,i2,j2)

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
 #   if !can_droop(bpd,i1,j1,i2,j2)
 #     return()
 #   end

 # assumes can_droop==true
    
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


function flat_droops(bpd)
   local n=size(bpd)[1]

   local dps = Vector{Matrix}([])

   for i1=1:n-1
     for j1=1:n-1
       for i2=i1+1:n
          for j2=j1+1:n
            if can_flat_droop(bpd,i1,j1,i2,j2)
              bpd2=flatten(droop(bpd,i1,j1,i2,j2))
              push!(dps,bpd2)
            end
          end
       end
     end
   end

   return(dps)
end


function sharp_droops(bpd)
   local n=size(bpd)[1]

   local dps = Vector{Matrix}([])

   for i1=1:n-1
     for j1=1:n-1
       for i2=i1+1:n
          for j2=j1+1:n
            if can_sharp_droop(bpd,i1,j1,i2,j2)
              bpd2=droop(bpd,i1,j1,i2,j2)
              push!(dps,bpd2)
            end
          end
       end
     end
   end

   return(dps)
end



# TO FIX: unsure about best lookup table to use here.
# Memory-wise, it seems bad to store all these values of a recursive function.
# But since they get met many times on different paths down, seems like there is a substantial speed-up.

#const bpd_list = Dict{Matrix, Vector}()

function all_below(bpd)

#   if haskey(bpd_list, bpd)
#     return bpd_list[bpd]
#   end

   local dps=all_droops(bpd)

   if length(dps)==0
#      bpd_list[bpd]=Vector{Matrix}([bpd])
#      return bpd_list[bpd]
     return [bpd]
   end

   local alldps = []

   for b in dps
     alldps=union(alldps,all_below(b))
   end

   prepend!(alldps,[bpd])

#   bpd_list[bpd]=alldps
#   return bpd_list[bpd]
   return alldps
end



@memoize function all_bpds(w)

  local bpd=Rothe(w)

  return all_below(bpd)

end



function is_sharp(bpd)
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

function is_flat(bpd)
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

 # assumes bpd is flat

   local dps=flat_droops(bpd)

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



@memoize function flat_bpds(w)

  local bpd=flatten(Rothe(w))

  return flat_below(bpd)

end




function sharps(w)

  filter(is_sharp, all_bpds(w))

end

function flats(w)

  filter(is_flat, all_bpds(w))

end


function drift_class(w,bpd)

  local bflat=flatten(bpd)

  dr = []
  
  for b in all_bpds(w)
    if flatten(b)==bflat
      push!(dr,b)
    end
  end

  return dr
end

function tableau_components(bpd)
# return labelled tableaux for a flat bpd

  local n=size(bpd)[1]

  if !is_flat(bpd)
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


#=

Using Nemo

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
  if !is_flat(bpd)
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

=#


# graphics
using Plots

function draw_se_elbow_curve(x1, y1, x3, y3)
    x2, y2 = x3, y1
    t = range(0, stop=1, length=100)
    x_vals = @. (1-t)^2 * x1 + 2*(1-t)*t * x2 + t^2 * x3
    y_vals = @. (1-t)^2 * y1 + 2*(1-t)*t * y2 + t^2 * y3
    return x_vals, y_vals
end

function draw_nw_elbow_curve(x1, y1, x3, y3)
    x2, y2 = x1, y3
    t = range(0, stop=1, length=100)
    x_vals = @. (1-t)^2 * x1 + 2*(1-t)*t * x2 + t^2 * x3
    y_vals = @. (1-t)^2 * y1 + 2*(1-t)*t * y2 + t^2 * y3
    return x_vals, y_vals
end


function draw_bpd(B,fn; img_size=(300,300))
    n, m = size(B)
    p = plot(; xlim=(0, n), ylim=(0, m), aspect_ratio=:equal, legend=false, grid=false, framestyle=:none, tick_direction=:none, size=img_size)
    
    for i = 1:n
        for j = 1:n
            y, x = n-i, j-1  # Transpose and invert the y-coordinate
            tile = B[i, j]
            if tile == "O"
                plot!([x, x+1, x+1, x, x], [y, y, y+1, y+1, y], linecolor=:orange, linewidth=0.5, seriestype=:shape, fillalpha=0)
            elseif tile == "+"
                plot!([x+0.5, x+0.5], [y, y+1], linecolor=:blue, linewidth=2)
                plot!([x, x+1], [y+0.5, y+0.5], linecolor=:blue, linewidth=2)
            elseif tile == "|"
                plot!([x+0.5, x+0.5], [y, y+1], linecolor=:blue, linewidth=2)
            elseif tile == "-"
                plot!([x, x+1], [y+0.5, y+0.5], linecolor=:blue, linewidth=2)
            elseif tile == "/"
                x_vals, y_vals = draw_se_elbow_curve(x+1, y+0.5, x+0.5, y)
                plot!(x_vals, y_vals, linecolor=:blue, linewidth=2)
            elseif tile == "%"
                x_vals, y_vals = draw_nw_elbow_curve(x+0.5, y+1, x, y+0.5)
                plot!(x_vals, y_vals, linecolor=:blue, linewidth=2)
            end
        end
    end
    plot!(; framestyle=:box, linecolor=:black, linewidth=3, ticks=nothing)
    savefig(fn)
end



function print_flat_bpds(w,fn,fmt)

  local ds=flat_bpds(w)
  local n=length(ds)

  for i=1:n
    draw_bpd(ds[i],string(fn,i,fmt))
  end

end




