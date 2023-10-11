# bpds
*julia code for bumpless pipe dreams*

This will eventually be bundled into a package for computing Schubert polynomials via drift polynomials and bumpless pipe dreams (BPDs).

For now, "bpds.jl" has functions which create all BPDs of a given permutation.  It can also create all the *flat* and *sharp* bpds, and perform droop and drop moves.

The app "bpd_app.jl" produces an interactive widget which displays BPDs and navigates by droop or drop moves among them, starting with the Rothe BPD.  To use it, make sure it is in the same directory as "bpds.jl", and that this directory is accessible to Julia.  Then open a Julia REPL and type
```include("bpd_app.jl")```

The "Rothe" button creates the Rothe BPD, and resets the window if other BPDs have been created.

The "Droops/Drops" toggle button determines which operation is triggered by clicking on a BPD.

Here are some screenshots.

<img src="https://github.com/pseudoeffective/bpds/assets/62109185/b4db5a87-b1b2-4c6c-ba16-c34730012840" height="400">
<img src="https://github.com/pseudoeffective/bpds/assets/62109185/1eaea7ba-7b47-44b3-9bb7-446d19667f08" height="400">
<img src="https://github.com/pseudoeffective/bpds/assets/62109185/6643d370-3313-4ec1-b0b4-3072ab90acb3" height="400">

Have fun!

NB: The app loads Gtk, which requires the Cairo package to be available to Julia (but not necessarily loaded).  If you get an error in the display referencing surface type, try running ```Pkg; Pkg.add("Cairo")``` and then restarting the REPL.
