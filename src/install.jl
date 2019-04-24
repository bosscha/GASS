## install package for Gaia.Clustering

using Pkg

Pkg.add("PyCall")
Pkg.add("DataFrames")
Pkg.add("Statistics")
Pkg.add("Distributions")
Pkg.add("Random")
Pkg.add("CSV")
Pkg.add("PyPlot")
Pkg.add("Query")
Pkg.add("LsqFit")
Pkg.add("FFTW")
Pkg.add("StatsBase")
Pkg.add("Distances")
Pkg.add("JLD")
Pkg.build("HDF5")


println("## Package installation for GASS done...")
