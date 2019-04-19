## Miscelaneous methods for GASS
##

## init population ...
function init_pop(obs , sub , wei , ga )
    arrcfg = CSV.read(obs.Array_Configuration_File, datarow=4 , header=["X" , "Y", "Z" , "diam" , "name"] , delim= " ")
    pop = population(arrcfg , obs , sub , wei , ga)

    return(pop)
end