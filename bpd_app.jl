# Widget for displaying Bumpless Pipe Dreams
# David Anderson, October 2023


using Gtk

include("./bpds.jl")

function parse_vector(input_str::String)
    # Remove brackets and split by comma
    elements = split(replace(input_str, "[" => "", "]" => ""), ",")
    # Convert to integers and return as a vector
    return parse.(Int, elements)
end


# Create the window, label, and entry
win = GtkWindow("BPD Maker", 400, 300)

label = Gtk.Label("Permutation w = ")
entry = Gtk.Entry()
GAccessor.text(entry, "[2,1,4,3]")

# Create toggle button
toggle_button = Gtk.ToggleButton("Drops")


# Create a button to generate the Rothe diagram
button = Gtk.Button("Rothe")

# Create a box for the images
action_grid = GtkGrid()

# Fill the grid with 4x5 event boxes

for i=1:5
  for j=1:4
    ebox = Gtk.EventBox()
    action_grid[i,j]=ebox
  end
end




# Pack the basic boxes

hbox = Gtk.Box(:h)
vbox = Gtk.Box(:v)

push!(hbox, label)
push!(hbox, entry)
push!(hbox, toggle_button)

push!(vbox, hbox)
push!(vbox, button)
push!(vbox, action_grid)



# Create a one-element array to hold the current bpd
current_bpd = Matrix{Any}(nothing,5,4)

# Create an array to keep track of the widgets added to the grid
grid_children = falses(5,4)

# Create an array to keep track of images
images = Matrix{Any}(undef,5,4)


# Tell the toggle button to display its state
signal_connect(toggle_button, "toggled") do widget
    is_flat = get_gtk_property(widget, :active, Bool)
    set_gtk_property!(widget, :label, is_flat ? "Drops" : "Droops")
end


# Connect the "Rothe" button to the action
signal_connect(button, "clicked") do widget

    for i=1:5
      for j=1:4
        if grid_children[i,j]
          delete!(action_grid[i,j],images[i,j])
          grid_children[i,j] = false
          current_bpd[i,j] = nothing
        end
      end
    end


    # Create an image widget
    image = Gtk.Image()

    push!(action_grid[1,1],image)


    user_input = get_gtk_property(entry, :text, String)
    w = parse_vector(user_input)
    bpd = Rothe(w)
    current_bpd[1,1]=bpd # Update the current bpd

    
    # Generate the image and update the Gtk.Image widget
    tmp_filename = tempname() * ".png"
    draw_bpd(bpd, tmp_filename)
    set_gtk_property!(image, :file, tmp_filename)
    rm(tmp_filename)

    grid_children[1,1] = true
    images[1,1] = image


    showall(action_grid)
    showall(vbox)

end




for i=1:5
  for j=1:4

# Connect a "button-press-event" to each Gtk.EventBox
signal_connect(action_grid[i,j], "button-press-event") do widget, event

    if current_bpd[i,j] != nothing  # Check if a bpd has been generated
        is_flat = get_gtk_property( toggle_button, :active, Bool )

        bpds = is_flat ? flat_drops(current_bpd[i,j]) : all_droops(current_bpd[i,j])
        num_images = minimum([length(bpds),20])
        
        # Compute the size dynamically based on the number of images
        computed_size = (400 ÷ sqrt(num_images), 400 ÷ sqrt(num_images))
        
        filenames = []
        for i=1:num_images
            tmp_filename = tempname() * "_$i.png"
            draw_bpd(bpds[i], tmp_filename; img_size=computed_size)
            push!(filenames, tmp_filename)
        end



    for i=1:5
      for j=1:4
        if grid_children[i,j]
          delete!(action_grid[i,j],images[i,j])
          grid_children[i,j]=false
          current_bpd[i,j]=nothing
        end
      end
    end


       k=1
       for i=1:5
         for j=1:4
           if k<=num_images
            img = Gtk.Image()
            fn=filenames[k]
            set_gtk_property!(img, :file, fn)
            push!(action_grid[i,j], img)
            images[i,j]=img
            grid_children[i,j]=true
            current_bpd[i,j]=bpds[k]
            rm(fn)  # Remove the temporary file
           end
           k +=1
          end
        end
       
        showall(action_grid)
    end
end

end

end

push!(win, vbox)

showall(win)
