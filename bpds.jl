function Rothe(w)
        local n=length(w)
        local r=Matrix{Union{Missing,String}}(missing,n,n)
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

   local dps = Set{Matrix}([])

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

   local dps = Set{Matrix}([])

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

   local dps = Set{Matrix}([])

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



const bpd_list = Dict{Matrix, Set}()

function all_below(bpd)

   if haskey(bpd_list, bpd)
     return bpd_list[bpd]
   end

   if length(all_droops(bpd))==0
      bpd_list[bpd]=Set{Matrix}([bpd])
      return bpd_list[bpd]
#     return Set{Matrix}([bpd])
   end

   local alldps = Set{Matrix}([])
   local dps=all_droops(bpd)
   for b in dps
     alldps=union(alldps,all_below(b))
   end

   push!(alldps,bpd)

   bpd_list[bpd]=alldps
   return bpd_list[bpd]
#   return alldps
end


function all_bpds(w)

  local bpd=Rothe(w)

  return( collect( all_below(bpd) ) )

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


const flat_bpd_list = Dict{Matrix, Set}()

function flat_below(bpd)

   if haskey(flat_bpd_list, bpd)
     return flat_bpd_list[bpd]
   end

   if length(flat_droops(bpd))==0
      flat_bpd_list[bpd]=Set{Matrix}([bpd])
      return flat_bpd_list[bpd]
#     return Set{Matrix}([bpd])
   end

   local alldps = Set{Matrix}([])
   local dps=flat_droops(bpd)
   for b in dps
     alldps=union(alldps,flat_below(b))
   end

   push!(alldps,bpd)

   flat_bpd_list[bpd]=alldps
   return flat_bpd_list[bpd]
#   return alldps
end

function flat_bpds(w)

  local bpd=Rothe(w)

  return( collect( flat_below(bpd) ) )

end




function sharps(w)

  filter(is_sharp, all_bpds(w))

end

function flats(w)

  filter(is_flat, all_bpds(w))

end



#=
# graphics
# using Plots

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


function draw_bpd(B,fn)
    n, m = size(B)
    p = plot(; xlim=(0, n), ylim=(0, n), aspect_ratio=:equal, legend=false, grid=false, framestyle=:none, tick_direction=:none)
    
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


# sample for producing a bunch of images
for i=1:240
       draw_bpd(tww3[i],string("./bpd_figs/w32187654bpd",i,".pdf"))
       end;


=#
