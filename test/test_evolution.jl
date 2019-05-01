## test to run a GASS optimization
##

rootdir = "/home/stephane/Science/ALMA/ArrayConfig/GASS"

push!(LOAD_PATH,"$rootdir/master/src")
using GASS

## directory
wdir    = "$rootdir/products"

cd(wdir)

function main(inpfile)
    println("# GASS subarray optimization")
    cfg = read_cfg(inpfile)
    
    res= gass_optimization(cfg)
    println(res)
end

pop= @time main("../master/data/GA_Parameters_O-3.txt.julia")
