## types for GASS
##

###########################
struct inputCfg
    file_parameters::String
    folder_results::String
    verbose::Bool
end

##########################
struct observation
    Array_Configuration_File::String
    Observatory_Latitude::Float64
    Source_Declination::Float64
    Source_Hour_Angle::Float64
    Subarray_Number::Int
end

############################
struct subarrayParameters
        Pads_Per_Subarray::Vector{Int}
        Subarray_Name::Vector{String}
        Subrange::Array{UnitRange{Int64},1}
        Spatial_Resolution::Vector{Float64}
        Maximum_Recoverable_Scale::Vector{Float64}
        Elongation::Vector{Float64}
        Sidelobe_Level::Vector{Float64}
end

###############################
struct weight
    Weight_Subarray::Vector{Float64}
    Weight_Spatial_Resolution::Vector{Float64}
    Weight_Maximum_Recoverable_Scale::Vector{Float64}
    Weight_Elongation::Vector{Float64}                  
    Weight_Sidelobe_Levels::Vector{Float64}
end

############################
struct GA
    Number_Iterations::Int
    Population_Size::Int
    Termination_Condition::Bool
    Threshold::Float64
    Mutation_Rate::Float64
    Tournament_Size::Int
    Number_Elitism::Int
end

############################
struct cfg
    arr::AbstractDataFrame 
    obs::observation
    sub::subarrayParameters
    wei::weight
    ga::GA
    inp::inputCfg
end

