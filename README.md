# bpds
*julia code for bumpless pipe dreams*

This will eventually be bundled into a package for computing Schubert polynomials via drift polynomials and bumpless pipe dreams (BPDs).  These are gadgets which enumerate terms in Schubert polynomials, developed by [Lam-Lee-Shimozono](https://arxiv.org/abs/1806.11233).

For now, "BpdBase.jl" has functions which create all BPDs of a given permutation.  It can also create all the *flat* and *sharp* bpds, and perform droop and drop moves.

The file "BpdDraw.jl" can produce nice-looking diagrams for BPDs, in several formats (PNG, PDF, etc.), depending on what backend the Plots package is using.

The app "BpdApp.jl" produces an interactive widget which displays BPDs and navigates by droop or drop moves among them, starting with the Rothe BPD.  To use it, make sure it is in the same directory as "BpdBase.jl" and "BpdDraw.jl", and that this directory is accessible to Julia.  Then open a Julia REPL and type
```include("BpdApp.jl")```

The "Rothe" button creates the Rothe BPD, and resets the window if other BPDs have been created.

The "Droops/Drops" toggle button determines which operation is triggered by clicking on a BPD.

Here are some screenshots.

<img src="https://github.com/pseudoeffective/bpds/blob/main/Rothe.jpg" height="400">
<img src="https://github.com/pseudoeffective/bpds/blob/main/Drops.jpg" height="400">
<img src="https://github.com/pseudoeffective/bpds/blob/main/Droops.jpg" height="400">

Have fun!

NB: The app loads Gtk, which requires the Cairo package to be available to Julia (but not necessarily loaded).  If you get an error in the display referencing surface type, try running ```Pkg; Pkg.add("Cairo")``` and then restarting the REPL.
