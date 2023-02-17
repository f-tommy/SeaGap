module SeaGap

# External packages
using Printf
using LinearAlgebra
using Random
using Dates
using Statistics
using Dierckx
using Optim
using Plots
using Base.Threads 
using Distributed
import DelimitedFiles
import PDFmerger
import GMT
import Distributions  

# Original functions
include("dateprocessing.jl")
include("ntdbasis.jl")
include("perturbation.jl")
include("read_gnssa.jl")
include("traveltime.jl")
include("anttena2tr.jl")
include("interpolate_gps.jl")
include("ll2xy.jl")
include("unixsort.jl")
include("obsdata_format.jl")
include("simple_inversion.jl")
include("LineFitting.jl")
include("running_filter.jl")
include("make_initial.jl")

# Original functions for positioning
include("ttres.jl")
include("denoise.jl")
include("kinematic_array.jl")
include("static_array.jl")
include("static_array_AICBIC.jl")
include("static_array_grad.jl")
include("static_array_TR.jl")
include("static_individual.jl")
include("static_array_mcmcgrad.jl")
include("static_array_mcmcgradc.jl")

# Original functions for postprocessing
include("plot_prof.jl")
include("plot_track.jl")
include("plot_kinematic_array.jl")
include("denoise_kinematic.jl")
include("position_kinematic.jl")
include("plot_ttres.jl")
include("plot_ntd.jl")
include("plot_AICBIC.jl")
include("plot_cormap.jl")
include("plot_histogram.jl")
include("plot_histogram2d.jl")
include("plot_histogram2d_each.jl")
include("plot_mcmcres.jl")
include("plot_mcmcparam.jl")
include("plot_mcmcparam_each.jl")
include("plot_gradmap.jl")
include("plot_ntdgrad.jl")

# Post-processing
include("convert_displacement.jl")
include("plot_displacement.jl")

# Test program
include("forward_test.jl")
end
