## to split the cfg files
## Format is :  "keyword : value"
function input_parameters(file)
    let
        param = []
        open(file,"r") do f
            for line in eachline(f)
                l = lstrip(line)
                if length(l) > 0 && l[1] != '#'
                    res = strip.(split(l,":"))
                    push!(param,res)
                end
            end
        end
        return(param)
    end
end

## parse input file
## 
function parse_input(input_param)
    let
        fileParam= "default.cfg" 
        saveFold=  "Results/"
        verbose= true
        
        for pair in input_param
            if pair[1] == "File_Parameters"
                fileParam = pair[2]
            elseif pair[1] == "Folder_Results"
                saveFold= pair[2]
            elseif pair[1] == "Display_Screen_Results" && pair[2] == "false"
                verbose= false
            end
        end
        res= inputCfg(fileParam , saveFold , verbose)
        return(res)    
    end
end

## parse parameter file for the subarray constraints
## No consistency check yet...
function parse_parameters(input_param)
    let
        obs= 0
        sub= 0
        wei= 0
        ga= 0
        
        Array_Configuration_File= "alma.cfg"
        Observatory_Latitude= -23.0262015
        Source_Declination= -30
        Source_Hour_Angle = -1 
        Subarray_Number= 4
        
        ## array of nsubarray elements
        Pads_Per_Subarray= []
        Weight_Subarray= []
        Subarray_Name= [] 
        Spatial_Resolution= []
        Maximum_Recoverable_Scale=
        Elongation= []
        Sidelobe_Level= []
        Weight_Spatial_Resolution= []
        Weight_Maximum_Recoverable_Scale= [] 
        Weight_Elongation= []                    
        Weight_Sidelobe_Levels= []
        
        ## GA parameters
        Number_Iterations = 100
        Population_Size= 100
        Termination_Condition= false
        Threshold= -0.05
        Mutation_Rate= 0.05
        Tournament_Size= 5
        Number_Elitism= 5
        
        for pair in input_param
            if pair[1] == "Array_Configuration_File"
                Array_Configuration_File= pair[2]
            elseif pair[1] == "Observatory_Latitude"
                Observatory_Latitude= parse(Float64, pair[2]) 
            elseif pair[1] == "Source_Declination"
                Source_Declination= parse(Float64, pair[2])
            elseif pair[1] == "Source_Hour_Angle"
                Source_Hour_Angle= parse(Float64, pair[2])
            elseif pair[1] == Subarray_Number
                Subarray_Number= parse(Int, pair[2])
            elseif pair[1] == "Number_Iterations"
                Number_Iterations= parse(Int, pair[2])
            elseif pair[1] == "Population_Size"
                Population_Size= parse(Int, pair[2])
            elseif pair[1] == "Termination_Condition"
                Termination_Condition= pair[2]=="true" ? true : false
            elseif pair[1] == "Threshold"
                Threshold=parse(Float64, pair[2]) 
            elseif pair[1] == "Mutation_Rate"
                Mutation_Rate=parse(Float64, pair[2])
            elseif pair[1] == "Tournament_Size"
                Tournament_Size= parse(Int, pair[2])
            elseif pair[1] == "Number_Elitism"
                Number_Elitism= parse(Int, pair[2])  
            elseif pair[1] == "Pads_Per_Subarray"
                Pads_Per_Subarray= map(x->(v = tryparse(Int,x); 
                        v==nothing ? 0.0 : v),split(pair[2],","))
            elseif pair[1] == "Weights_Subarray"
                  Weight_Subarray= map(x->(v = tryparse(Float64,x);
                        v==nothing ? 0.0 : v),split(pair[2],","))  
            elseif pair[1] == "Subarray_Name" 
                Subarray_Name= split(pair[2],",")
            elseif pair[1] == "Spatial_Resolution"
                Spatial_Resolution= map(x->(v = tryparse(Float64,x); 
                        v==nothing ? 0.0 : v),split(pair[2],","))
            elseif pair[1] == "Maximum_Recoverable_Scale"
                Maximum_Recoverable_Scale= map(x->(v = tryparse(Float64,x); 
                        v==nothing ? 0.0 : v),split(pair[2],","))       
            elseif pair[1] == "Elongation"
                Elongation= map(x->(v = tryparse(Float64,x);
                        v==nothing ? 0.0 : v),split(pair[2],","))
            elseif pair[1] == "Sidelobe_Level"
               Sidelobe_Level= map(x->(v = tryparse(Float64,x);
                        v==nothing ? 0.0 : v),split(pair[2],",")) 
            elseif pair[1] == "Weight_Spatial_Resolution"
                Weight_Spatial_Resolution= map(x->(v = tryparse(Float64,x);
                        v==nothing ? 0.0 : v),split(pair[2],","))
            elseif pair[1] == "Weight_Maximum_Recoverable_Scale" 
                Weight_Maximum_Recoverable_Scale= map(x->(v = tryparse(Float64,x);
                        v==nothing ? 0.0 : v),split(pair[2],","))
            elseif pair[1] == "Weight_Elongation" 
                Weight_Elongation= map(x->(v = tryparse(Float64,x);
                        v==nothing ? 0.0 : v),split(pair[2],","))
            elseif pair[1] == "Weight_Sidelobe_Levels"
                Weight_Sidelobe_Levels=map(x->(v = tryparse(Float64,x);
                        v==nothing ? 0.0 : v),split(pair[2],","))
            end
        end
        ### setting the struct.
        obs= observation(Array_Configuration_File, Observatory_Latitude , Source_Declination,Source_Hour_Angle,
            Subarray_Number)
        sub= subarrayParameters(Pads_Per_Subarray, Subarray_Name , Spatial_Resolution,
            Maximum_Recoverable_Scale , Elongation, Sidelobe_Level)
        wei= weight(Weight_Subarray,Weight_Spatial_Resolution,Weight_Maximum_Recoverable_Scale, 
            Weight_Elongation,Weight_Sidelobe_Levels)
        ga= GA(Number_Iterations , Population_Size,Termination_Condition ,Threshold ,
            Mutation_Rate , Tournament_Size ,Number_Elitism)
        
   return(obs, sub , wei , ga)     
end
end


## Read all the inputs from the input file
function read_input_cfg(inpfile)
    res= input_parameters(inpfile)
    inpcfg= parse_input(res)
    
    ## parameters inputs
    res= input_parameters(inpcfg.file_parameters)
    paramcfg= parse_parameters(res)
    
    return(inpcfg , paramcfg)
end
