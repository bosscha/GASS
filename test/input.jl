## test the initialisation of the optimisation

rootdir = "/home/stephane/Science/ALMA/ArrayConfig/GASS"

push!(LOAD_PATH,"$rootdir/master/src")
Using GASS

## directory
datadir = "$rootdir/master/data"
wdir    = "$rootdir/products"
plotdir = "$rootdir/products/test"

cd(wdir)

macro main(inpfile)
    res= input_parameters(inpfile)
    inpcfg= parse_input(res)
    println(inpcfg)
    
    ## parameters inputs
    res= input_parameters(inpcfg[:file_parameters])
    paramcfg= parse_param(res)
    println(paramcfg)
end

@main("../master/data/GA_Inputs_O-10.txt.julia")