## full test to run a GASS optimization with the output
##

rootdir = "/home/stephane/Science/ALMA/ArrayConfig/GASS"

push!(LOAD_PATH,"$rootdir/master/src")
using GASS

## directory
wdir    = "$rootdir/products"

cd(wdir)

function main(inpfile)
    cfg = read_cfg(inpfile , verbose=true)
    
    res= full_run("test_" , cfg , true)
end

pop= @time main("../master/data/GA_Inputs_O-3.txt.julia")