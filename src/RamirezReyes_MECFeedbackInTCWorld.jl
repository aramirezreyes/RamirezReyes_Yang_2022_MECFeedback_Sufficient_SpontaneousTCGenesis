module RamirezReyes_MECFeedbackInTCWorld

using Distributed: @everywhere, @spawnat, fetch
using Unitful: @u_str, unit, ustrip
using ImageSegmentation
using AvailablePotentialEnergyFramework
using JLD
using NCDatasets,NetCDF
using Statistics

include("APE_calculation_from_sam_output.jl")
include("Composite_creation_from_sam_output_and_diagnostics_in_chunks_9hpa.jl")
include("count_updrafts_and_inflow.jl")
export getapeanalysis,getapeanalysis_nosmoothing,get_composites, getapeanalysis_last15days, get_composites_in_chunks
export readfile_and_count_updrafts_and_downdrafts, read_file_and_compute_updraft_ratio_and_mean_inflow

end
