## methods to output the GASS results, mainly in CASA format
##

### save an array in CASA format using cfg mapping
###
function save_CASAarr(fileprefix::String, cfg::cfg, bestsub)    
    headerCasa= "# observatory=ALMA\n# coordsys=LOC (local tangent plane)\n# x y z diam pad#\n"
    
    for i in 1:length(bestsub)
        arrStr= string(headerCasa)
        subname= cfg.sub.Subarray_Name[i]

        for j in 1:length(bestsub[i])
            arrStr *= @sprintf("%f %f %f %2.1f %s \n",cfg.arr.X[bestsub[i][j]], cfg.arr.Y[bestsub[i][j]] , 
                cfg.arr.Z[bestsub[i][j]] , cfg.arr.diam[bestsub[i][j]] , cfg.arr.name[bestsub[i][j]])
        end

        filename= fileprefix*"_"*subname*".cfg"đ
        f= open(filename, "w")
        print(f,arrStr)
        close(f)
        println("## $filename written...")
    end
end

## sort population wrt. score
##
function sort_population(p::population)
    s= sortperm(p.score, rev=true)
    psorted= population(p.age, p.subarr[s,:], p.fitness[s,:], p.score[s], p.param[s,:])
    return(psorted)
end

## full run with JLD output ..
##
function full_run(fileprefix::String, cfg::cfg, jld= false)
    res= gass_optimization(cfg)

    lastPopulation= cfg.ga.Number_Iterations
    psort= sort_population(res[lastPopulation])
    save_CASAarr(fileprefix, cfg, psort.subarr[1,:])
    
    if jld
       JLD.save(fileprefix*".jld", "gass", res)
       println("## $fileprefix.jld written..")
    end
    
    return(res)
end
