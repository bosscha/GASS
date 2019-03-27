## Use the astropy python package to read a votable.

module GASS

using PyCall
using DataFrames
using Printf
using LsqFit , FFTW
using Statistics , Distributions ,Random
using StatsBase

import CSV
import PyPlot

include("types.jl")
export inputCfg , observation , subarrayParameters , weight , GA , cfg

include("inputcfg.jl")
export input_parameters , parse_input , parse_parameters , read_cfg

include("utils.jl")
export init_pop

include("simobserve.jl")
export calc_baselines , calc_uv , calc_dirtybeam , fit_beam , calc_mrs

include("evolution.jl")
export fitness_subarray , create_population , get_elitism , get_parents , get_crossover1 , 
  get_mutation , get_evolution , gass_optimization

end