## Genetic Algorithm 
##


## fitness function for a subarray
## cfg: GASS parameters
## subarrid: subarray id (int)
## subind: subarray indices in the main Arr.
##
function fitness_subarray(cfg, subarrid, subind)
    subarr = cfg.arr[subind,:]
    
    bl= calc_baselines(subarr)
    uv= calc_uv(bl, cfg.obs.Source_Hour_Angle ,  cfg.obs.Source_Declination)
    h , dr=  calc_dirtybeam(uv , 255, 127, robust=0.5)
    b= fit_beam(h , dr)
    mrs= calc_mrs(uv)
    
    #@printf("## subarray fitness \n")
    #@printf("## beam:")
    #println(b)
    #@printf("## MRS: %3.3f \n", mrs)
    
    res= 0
    res += cfg.wei.Weight_Spatial_Resolution[subarrid]*abs(b.ar-cfg.sub.Spatial_Resolution[subarrid])
    res += cfg.wei.Weight_Elongation[subarrid]*abs(b.e-cfg.sub.Elongation[subarrid])*sign(b.e-
        cfg.sub.Elongation[subarrid])
    res += cfg.wei.Weight_Sidelobe_Levels[subarrid]*abs(b.sidelobe-
        cfg.sub.Sidelobe_Level[subarrid])*sign(b.sidelobe-cfg.sub.Sidelobe_Level[subarrid])
    res += cfg.wei.Weight_Maximum_Recoverable_Scale[subarrid]*abs(mrs-
        cfg.sub.Maximum_Recoverable_Scale[subarrid])*sign(cfg.sub.Maximum_Recoverable_Scale[subarrid]-mrs)
    
    return(res , [b , mrs])
end


## Creation of a population
## 
function create_population(cfg)
    
    pop= Array{Array{Int,1},2}(undef,cfg.ga.Population_Size,  cfg.obs.Subarray_Number)
    fitness= Array{Float64,2}(undef,cfg.ga.Population_Size,  cfg.obs.Subarray_Number)
    score=  Array{Float64,1}(undef,cfg.ga.Population_Size)
    paramsub= Array{Dict{String,Float64},2}(undef,cfg.ga.Population_Size,  cfg.obs.Subarray_Number)
    
    subind= collect(1:cfg.obs.Antenna_Number)
    
    for i in 1:cfg.ga.Population_Size
        Random.shuffle!(subind)
        for j in 1:cfg.obs.Subarray_Number
            pop[i,j]= subind[cfg.sub.Subrange[j]]
            fitness[i,j] , res= _fitness_subarray(cfg, j, pop[i,j])
            println(i," ",j," ",fitness[i,j])
            println(res)
            paramsub[i,j]= Dict("ar"=>res[1].ar,"e"=>res[1].e, "sidelobe"=>res[1].sidelobe, "mrs"=>res[2])
        end
        score[i]= -sum(cfg.wei.Weight_Subarray[:] .* fitness[i,:])
    end
    
    pop0= population(0, pop, fitness, score , paramsub)
    
    return(pop0)
end

### get best parents...
###
function get_elitism(cfg, pop::population)
    isort= reverse(sortperm(pop.score))
    
    elit= []
    for i in 1:cfg.ga.Number_Elitism
        push!(elit, [pop.subarr[isort[i],:], pop.fitness[isort[i],:], pop.score[isort[i]]])
    end
    return(elit)
end

## ... using tournament selection
function get_parents(cfg, pop::population)
    tour= cfg.ga.Tournament_Size
    
    popran= Random.shuffle!(collect(1:cfg.ga.Population_Size))
    isort= reverse(sortperm(pop.score[popran[1:tour]]))
    iparent1= popran[isort[1]]
  
    iparent2= iparent1
    while iparent1==iparent2
        popran= Random.shuffle!(collect(1:cfg.ga.Population_Size))
        isort= reverse(sortperm(pop.score[popran[1:tour]]))
        iparent2= popran[isort[1]] 
    end 
    return(iparent1 , iparent2)
end

#### auxiliary function
function flatten_subarray(cfg, subarr)
    subflat= zeros(Int, cfg.obs.Antenna_Number)
    
    counter= 1
    for i in 1:cfg.obs.Subarray_Number
        for j in 1:cfg.sub.Pads_Per_Subarray[i]
            subflat[counter]= subarr[i][j]
            counter += 1
        end
    end
    return(subflat)
end

### recreate the subarray array..
function wrap_subarray(cfg, sub)
    subarray= Array{Array{Int,1}}(undef,0)

    for i in 1:cfg.obs.Subarray_Number
        push!(subarray, sub[cfg.sub.Subrange[i]])
    end
    
    return(subarray)
end

### crossover (#1) using the two parents
###
function get_crossover1(cfg, pop::population , parents)
    nmax= cfg.obs.Antenna_Number
    p1= _flatten_subarray(cfg, pop.subarr[parents[1],:])
    p2= _flatten_subarray(cfg, pop.subarr[parents[2],:])
    child= zeros(Int, nmax)

    pivot= rand(1:nmax,2)
    pivot1= minimum(pivot) ; pivot2= maximum(pivot)
    child[pivot1:pivot2] .= p1[pivot1:pivot2]

    irange1= 1:pivot1-1 ; irange2= pivot2+1:nmax
    
    cnterp2= 1
    for i in irange1
        while (p2[cnterp2] in child)
            cnterp2 += 1
        end
        child[i]= p2[cnterp2]
    end
    for i in irange2
        while (p2[cnterp2] in child)
            cnterp2 += 1
        end
        child[i]= p2[cnterp2]
    end  
    
    # println(child)
    c= _wrap_subarray(cfg, child)
    # println(c)
    return(child)
end

### mutation
### loop over all alleles and swap two if p is true
### 
function get_mutation(cfg, child)
    nmax= cfg.obs.Antenna_Number
    si= copy(child)
    
    pmutation= 1-cfg.ga.Mutation_Rate
    rng = Random.MersenneTwister()
    rd= Random.rand(rng, nmax)
    # println(rd)
    for i in 1:nmax
        if rd[i] > pmutation
            iswap= Random.rand(rng, big.(1:nmax))
            value= si[iswap]
            si[iswap]= si[i]
            si[i]= value
        end
    end
    return(si)
end

###
### Evolve a population
###
function get_evolution(cfg, pi::population)
    npop= cfg.ga.Population_Size
    age= pi.age + 1
    nelit= cfg.ga.Number_Elitism
    
    ## Elitism
    pelit=  get_elitism(cfg, pi)
    
    ## parent selection, crossover and mutation
    ncross= npop-nelit
    
    crosspop= []
    for i in 1:ncross
        parent= get_parents(cfg, pi)
        child=  get_crossover1(cfg, pi , parent)
        mutated= get_mutation(cfg, child)
        push!(crosspop, mutated)
        println(parent)
        println(child)
    end

    ## create an evolved population
    pinew= Array{Array{Int,1},2}(undef,cfg.ga.Population_Size,  cfg.obs.Subarray_Number)
    fitness= Array{Float64,2}(undef,cfg.ga.Population_Size,  cfg.obs.Subarray_Number)
    score=  Array{Float64,1}(undef,cfg.ga.Population_Size)
    paramsub= Array{Dict{String,Float64},2}(undef,cfg.ga.Population_Size,  cfg.obs.Subarray_Number)
    
    for i in 1:nelit
        for j in 1:cfg.obs.Subarray_Number
            pinew[i,j]= pelit[i][1][j]
            fitness[i,j], res= _fitness_subarray(cfg, j, pinew[i,j])
            println(i," ",j," ",fitness[i,j])
        end
        score[i]= -sum(cfg.wei.Weight_Subarray[:] .* fitness[i,:])
    end    
    
    println("fitness of crossover..")
    for i in nelit+1:npop
        println(crosspop[i-nelit])
        pwrap= _wrap_subarray(cfg, crosspop[i-nelit])
        println(pwrap)
        for j in 1:cfg.obs.Subarray_Number
            pinew[i,j]= pwrap[j]
            fitness[i,j] , res= _fitness_subarray(cfg, j, pinew[i,j])
            paramsub[i,j]= Dict("ar"=>res[1].ar,"e"=>res[1].e, "sidelobe"=>res[1].sidelobe, "mrs"=>res[2])
            println(i," ",j," ",fitness[i,j])
        end
        score[i]= -sum(cfg.wei.Weight_Subarray[:] .* fitness[i,:])
    end
    
    popnew= _population(age, pinew, fitness, score, paramsub)
    
    return(popnew)
end

