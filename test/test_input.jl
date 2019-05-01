## test the initialisation of the optimisation

rootdir = "/home/stephane/Science/ALMA/ArrayConfig/GASS"

push!(LOAD_PATH,"$rootdir/master/src")
using GASS

## directory
wdir    = "$rootdir/products"

cd(wdir)

function main(inpfile)
    res= read_cfg(inpfile)
    println(res)
end

main("../master/data/GA_Parameters_O-3.txt.julia")