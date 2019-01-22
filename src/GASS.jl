## Use the astropy python package to read a votable.

module GASS

using PyCall
using DataFrames
using Printf
import CSV

using Statistics , Distributions ,Random

import PyPlot

include("types.jl")
export inputCfg , observation , subarrayParameters , weight , GA , cfg

include("inputcfg.jl")
export input_parameters , parse_input , parse_parameters , read_cfg

include("utils.jl")
export init_pop

end