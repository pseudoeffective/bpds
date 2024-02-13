# BPD posets
# David Anderson
# February 2024

include("BpdBase.jl")


function droop_poset(w)
# form list of pairs (i,j) where i-th BPD droops to j-th one

  cover_rels = []

  bpds=collect(all_bpds(w))

  bpds=map( b -> b.m, bpds )

  for i=1:length(bpds)
    drps = all_droops( BPD(bpds[i]) )
    drps = map( b -> b.m, drps )
    for b in drps
      j=findfirst(==(b),bpds)
      push!(cover_rels, (i,j))
    end
  end

  return cover_rels

end


function drift_poset(w)
# form list of pairs (i,j) where i-th BPD drifts to j-th one

  cover_rels = []

  bpds=collect(all_bpds(w))

  bpds=map( b -> b.m, bpds )

  for i=1:length(bpds)
    drps = all_drips( BPD(bpds[i]) )
    drps = map( b -> b.m, drps )
    for b in drps
      j=findfirst(==(b),bpds)
      push!(cover_rels, (i,j))
    end
  end

  return cover_rels

end


function exp_poset(po, fn)
    open(fn, "w") do file
#        for pair in po
#            print(file, pair)
#        end
      println(file, "[", join(po, ", "), "]")
    end
end
