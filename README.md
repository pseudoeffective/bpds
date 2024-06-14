# bpds
*julia widget for bumpless pipe dreams*


The file "BpdWidget.jl" produces an interactive widget which displays BPDs and navigates by droop or drop moves among them, starting with the Rothe BPD.  To use it, first load the [SchubertPolynomials](https://github.com/pseudoeffective/SchubertPolynomials.jl) package.  Make sure the file "BpdWidget.jl" is in a directory is accessible to Julia.  Then type
```julia-repl
julia> include("BpdApp.jl")
```

The "Rothe" button creates the Rothe BPD, and resets the window if other BPDs have been created.

The "Droops/Drops" toggle button determines which operation is triggered by clicking on a BPD.

Here are some screenshots.

<img src="https://github.com/pseudoeffective/bpds/blob/main/Rothe.jpg" height="400">
<img src="https://github.com/pseudoeffective/bpds/blob/main/Drops.jpg" height="400">
<img src="https://github.com/pseudoeffective/bpds/blob/main/Droops.jpg" height="400">

Have fun!

NB: The app loads Gtk, which requires the Cairo package to be available to Julia (but not necessarily loaded).  If you get an error in the display referencing surface type, try running ```Pkg; Pkg.add("Cairo")``` and then restarting the REPL.

