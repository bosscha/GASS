## Use the astropy python package to read a votable.

module GASS

using PyCall
using DataFrames
using Printf
import CSV
using LsqFit , FFTW

using Statistics , Distributions ,Random
using StatsBase

import PyPlot

include("types.jl")
export inputCfg , observation , subarrayParameters , weight , GA , cfg

include("inputcfg.jl")
export input_parameters , parse_input , parse_parameters , read_cfg

include("utils.jl")
export init_pop

include("simobserve.jl")
export calc_baselines , calc_uv , calc_dirtybeam , fit_beam

end