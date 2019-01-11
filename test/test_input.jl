## test the initialisation of the optimisation

rootdir = "/home/stephane/Science/ALMA/ArrayConfig/GASS"

push!(LOAD_PATH,"$rootdir/master/src")
using GASS

## directory
wdir    = "$rootdir/products"

cd(wdir)

macro main(inpfile)
    res= read_cfg(inpfile)
    println(res)
end

@main("../master/data/GA_Inputs_O-10.txt.julia")