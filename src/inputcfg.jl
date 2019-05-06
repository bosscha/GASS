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
        arraycfg= 0 
        obs= 0
        sub= 0
        wei= 0
        ga= 0
        
        Array_Configuration_File= "alma.cfg"
        Result_Folder= "."
        Display_Verbose= true
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
            elseif pair[1] == "Result_Folder"
                Result_Folder= pair[2]
            elseif pair[1] == "Display_Verbose"
                Display_Verbose=  pair[2]=="true" ? true : false
            elseif pair[1] == "Observatory_Latitude"
                Observatory_Latitude= parse(Float64, pair[2])    
            elseif pair[1] == "Source_Declination"
                Source_Declination= parse(Float64, pair[2])
            elseif pair[1] == "Source_Hour_Angle"
                Source_Hour_Angle= parse(Float64, pair[2])
            elseif pair[1] == "Subarray_Number"
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
                for i in 1:length(Subarray_Name)
                  Subarray_Name[i]= strip(Subarray_Name[i])
                end
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
        
        arrcfg = CSV.read(convert(String,Array_Configuration_File), datarow=4 , header=["X" , "Y", "Z" , "diam" , "name"] ,  ignorerepeated=true, delim= " ")
           
        subrange = []
        start= 1
            for npad in Pads_Per_Subarray
                push!(subrange,start:start+npad-1)
                start= start+npad
            end
        
        Antenna_Number= size(arrcfg,1)
        
        ### setting the struct.
        obs= observation(Array_Configuration_File, Result_Folder, Display_Verbose, Observatory_Latitude , Source_Declination,Source_Hour_Angle,
            Antenna_Number ,Subarray_Number)
        sub= subarrayParameters(Pads_Per_Subarray, Subarray_Name , subrange, Spatial_Resolution,
            Maximum_Recoverable_Scale , Elongation, Sidelobe_Level)
        wei= weight(Weight_Subarray,Weight_Spatial_Resolution,Weight_Maximum_Recoverable_Scale, 
            Weight_Elongation,Weight_Sidelobe_Levels)
        ga= GA(Number_Iterations , Population_Size,Termination_Condition ,Threshold ,
            Mutation_Rate , Tournament_Size ,Number_Elitism)
          
        ### population setting
        res= cfg(arrcfg , obs , sub , wei , ga)
        
   return(res)     
end
end

## Normalize the weights
##
function normalize_weight!(cfg)
   k1= sum(cfg.wei.Weight_Spatial_Resolution) + sum(cfg.wei.Weight_Maximum_Recoverable_Scale) + 
     sum(cfg.wei.Weight_Elongation) + sum(cfg.wei.Weight_Sidelobe_Levels)
   cfg.wei.Weight_Spatial_Resolution= cfg.wei.Weight_Spatial_Resolution ./ k1
   cfg.wei.Weight_Maximum_Recoverable_Scale= cfg.wei.Weight_Maximum_Recoverable_Scale ./ k1
   cfg.wei.Weight_Elongation= cfg.wei.Weight_Elongation ./ k1
   cfg.wei.Weight_Sidelobe_Levels= cfg.wei.Weight_Sidelobe_Levels ./ k1
  
  k2= sum(cfg.wei.Weight_Subarray)
  cfg.wei.Weight_Subarray= cfg.wei.Weight_Subarray ./ k2
end

## main function to read the cfg and check it.
function read_cfg(inpfile, check=true)
    res= input_parameters(inpfile)
    cfg= parse_parameters(res)
    normalize_weight!(cfg)
    
    if cfg.obs.Display_Verbose
      @printf("## Input Parameters for GASS \n")
      @printf("### Configuration file: %s \n", cfg.obs.Array_Configuration_File)
      @printf("### Result folder: %s \n", cfg.obs.Result_Folder)
      @printf("### Obs. Latitude: %3.3f \n", cfg.obs.Observatory_Latitude)
      @printf("### Source Declination: %3.1f \n", cfg.obs.Source_Declination)
      @printf("### HA: %3.1f \n", cfg.obs.Source_Hour_Angle)
      @printf("### Antenna number: %d \n", cfg.obs.Antenna_Number)
      @printf("### Subarray number: %d \n", cfg.obs.Subarray_Number)
      @printf("##\n")
      
      @printf("## Subarray Parameters\n")
      println("### Pads per subarray: ",cfg.sub.Pads_Per_Subarray)
      println("### Name: ",cfg.sub.Subarray_Name)
      println("### AR: ",cfg.sub.Spatial_Resolution)
      println("### MRS: ",cfg.sub.Maximum_Recoverable_Scale)
      println("### elongation: ",cfg.sub.Elongation)
      println("### sidelobe: ",cfg.sub.Sidelobe_Level)
      @printf("##\n")
      
      @printf("## GA parameters\n")
      @printf("### Iterations: %d \n",cfg.ga.Number_Iterations)
      @printf("### Population size: %d \n",cfg.ga.Population_Size)
      @printf("### Mutation rate: %3.3f \n",cfg.ga.Mutation_Rate)
      @printf("### Tournament size: %d \n",cfg.ga.Tournament_Size)
      @printf("### Elitism: %d \n",cfg.ga.Number_Elitism)
      @printf("##\n")
      
      @printf("## Weights\n")
      println("### Subarray weights: ",cfg.wei.Weight_Subarray)
      println("### AR weights: ",cfg.wei.Weight_Spatial_Resolution)
      println("### MRS weights: ",cfg.wei.Weight_Maximum_Recoverable_Scale)
      println("### elongation weights: ",cfg.wei.Weight_Elongation)
      println("### sidelobe weights: ",cfg.wei.Weight_Sidelobe_Levels)
      @printf("##\n")
    end
    
    if check
      check_consistency(cfg)
    end
    
    return(cfg)
end

function check_consistency(cfg::cfg)
    nsubpads= sum(cfg.sub.Pads_Per_Subarray)
    npads= nrow(cfg.arr)
    if nsubpads != npads
        error("##Error: Subarray pads are not equal to the total of pads.")
    end
    return(true)
end

