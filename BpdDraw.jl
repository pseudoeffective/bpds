# Tools for drawing BPDs in Julia
# David Anderson, October 2023.

using Plots

include("./BpdBase.jl")


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


