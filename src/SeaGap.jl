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
include("pos_array_each.jl")
include("pos_array_all.jl")
include("pos_array_all_AICBIC.jl")
include("pos_array_pvg.jl")
include("pos_array_TR.jl")
include("pos_single.jl")
include("pos_array_mcmcpvg.jl")
include("pos_array_mcmcpvgc.jl")

# Original functions for postprocessing
include("plot_prof.jl")
include("plot_track.jl")
include("plot_array_each.jl")
include("denoise_each.jl")
include("position_each.jl")
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
