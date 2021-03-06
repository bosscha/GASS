## full test to run a GASS optimization with the output
##

rootdir = "/home/stephane/Science/ALMA/ArrayConfig/GASS"

push!(LOAD_PATH,"$rootdir/master/src")
using GASS

## directory
wdir    = "$rootdir/products"

cd(wdir)

function main(inpfile)
    println("# GASS subarray optimization")
    println("## Parameter file: $inpfile")
    cfg = read_cfg(inpfile)
    
    res= full_run("testing" , cfg , true)
end

pop= @time main("../master/data/GA_Parameters_testing.txt.julia")
