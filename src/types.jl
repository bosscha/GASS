## types for GASS
##

###########################
mutable struct  inputCfg
    file_parameters::String
    folder_results::String
    verbose::Bool
end

##########################
mutable struct  observation
    Array_Configuration_File::String
    Observatory_Latitude::Float64
    Source_Declination::Float64
    Source_Hour_Angle::Float64
    Antenna_Number::Int
    Subarray_Number::Int
end

############################
mutable struct  subarrayParameters
        Pads_Per_Subarray::Vector{Int}
        Subarray_Name::Vector{String}
        Subrange::Array{UnitRange{Int64},1}
        Spatial_Resolution::Vector{Float64}
        Maximum_Recoverable_Scale::Vector{Float64}
        Elongation::Vector{Float64}
        Sidelobe_Level::Vector{Float64}
end

###############################
mutable struct  weight
    Weight_Subarray::Vector{Float64}
    Weight_Spatial_Resolution::Vector{Float64}
    Weight_Maximum_Recoverable_Scale::Vector{Float64}
    Weight_Elongation::Vector{Float64}                  
    Weight_Sidelobe_Levels::Vector{Float64}
end

############################
mutable struct GA
    Number_Iterations::Int
    Population_Size::Int
    Termination_Condition::Bool
    Threshold::Float64
    Mutation_Rate::Float64
    Tournament_Size::Int
    Number_Elitism::Int
end

############################
############################
mutable struct cfg
    arr::AbstractDataFrame 
    obs::observation
    sub::subarrayParameters
    wei::weight
    ga::GA
    inp::inputCfg
end

############################
## size are in arcsecs
struct synthbeam
    bx::Float64       ## bx inarcsec
    by::Float64       ## by i arcsec
    ar::Float64       ## angular resolution in arcsec
    e::Float64        ## excentricity
    sidelobe::Float64 ## sidelobe levels in percentage
end

#############################
#############################
mutable struct population
    age::Int32                        ## age
    subarr::Array{Array{Int,1},2}     ## {Population,Subindices}
    fitness::Array{Float64,2}         ## fitness of each subarray
    score::Array{Float64,1}           ## global score of each set
    param::Array{Dict{String,Float64},2}    ## beam and mrs for each subarray    
end


