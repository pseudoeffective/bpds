# bpds
*julia code for bumpless pipe dreams*

This will eventually be bundled into a package for computing Schubert polynomials via drift polynomials and bumpless pipe dreams (BPDs).  These are gadgets which enumerate terms in Schubert polynomials, developed by [Lam-Lee-Shimozono](https://arxiv.org/abs/1806.11233).

The file "BpdBase.jl" has functions which create all BPDs of a given permutation.  It can also create all the *flat* and *sharp* BPDs, perform droop and drop moves, and do the analogous things for the K-theory versions.

The file "BpdSchub.jl" has several methods for producing (double) Schubert polynomials.  For example,
```julia
julia> include("BpdSchub.jl");

julia> w = [1,4,3,2];

# set the ambient ring, this one just in x variables
julia> Rx = xy_ring(4,0)[1];

julia> sp = schub_poly(w,Rx)
x1*x2 + x1*x3 + x2*x3

# for the double Schubert polynomial, change the ring
julia> Rxy = xy_ring(4,4)[1];

julia> sp = schub_poly(w,Rxy)
x1*x2 + x1*x3 + x1*y1 + x1*y2 + x2*x3 + x2*y1 + x2*y2 + x3*y1 + x3*y2 + y1^2 + y1*y2 + y2^2
```

There is a function for quickly counting the terms of a Schubert polynomial:
```julia
julia> w=[1,4,3,2,10,9,8,7,6,5];

julia> nschub(w)
4424420
```
For more info, see the help files:
```julia
julia> ?schub_poly

julia> ?groth_poly
```

The file "SSYT.jl" has methods for generating Young tableaux, as well as (double, flagged) Schur polynomials.  It is loaded with "BpdSchub.jl".
```julia
julia> la = [2,1]; ff=2;

julia> ssyt( la, ff )
2-element Vector{Tableau}:
 
1 1 
2 

 
1 2 
2 

julia> ff=[2,3];

julia> R = xy_ring(5,0)[1];

julia> sp = schur_poly( la, ff, R )
x1^2*x2+x1*x2^2

# the argument ff can be a vector of integers, for a flagged Schur polynomial:
julia> sp = schur_poly( la, [2,3], R)
x1^2*x2 + x1^2*x3 + x1*x2^2 + x1*x2*x3 + x2^2*x3
```

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

*Thanks to [Anders Buch](https://sites.math.rutgers.edu/~asbuch/) for convincing me to try out Julia.*
