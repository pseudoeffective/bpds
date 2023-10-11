# bpds
julia code for bumpless pipe dreams

This will eventually be bundled into a package for computing Schubert polynomials via drift polynomials and bumpless pipe dreams (BPDs).

For now, "bpds.jl" has functions which create all BPDs of a given permutation w.  It can also create all the *flat* and *sharp* bpds, and perform droop and drop moves.

The app "bpd_app.jl" produces an interactive widget which displays BPDs and navigates by droop or drop moves among them, starting with the Rothe BPD.  To use it, open a Julia REPL and type
```include("bpd_app.jl")```

Here are some screenshots.

The app loads Gtk, which requires the Cairo package to be available to Julia (but not necessarily loaded).  If you get an error in the display referencing surface type, try running ```Pkg; Pkg.add("Cairo")``` and then restarting the REPL.
