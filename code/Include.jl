# =============================================================================
# Include.jl
# Centralized environment setup, package imports, and source loading for the
# hybrid SIM composition experiments.
# =============================================================================

const _ROOT = @__DIR__
const _PATH_TO_SRC      = joinpath(_ROOT, "src")
const _PATH_TO_DATA     = joinpath(_ROOT, "data")
const _PATH_TO_FIGS     = joinpath(_ROOT, "figs")
const _PATH_TO_SCRIPTS  = joinpath(_ROOT, "scripts")
const _PATH_TO_CONFIG   = joinpath(_ROOT, "config.toml")

using Pkg
Pkg.activate(_ROOT)
if !isfile(joinpath(_ROOT, "Manifest.toml"))
    # Order matters: VLQuantitativeFinancePackage depends on JumpHMM, so JumpHMM
    # must be registered in this environment first.
    Pkg.add(url = "https://github.com/varnerlab/JumpHMM.jl.git")
    Pkg.add(url = "https://github.com/varnerlab/VLQuantitativeFinancePackage.jl.git")
    Pkg.resolve()
    Pkg.instantiate()
end

using VLQuantitativeFinancePackage
import JumpHMM
using JumpHMM: JumpHiddenMarkovModel, SingleIndexModel, HybridSingleIndexModel,
               StudentTCopula, GaussianCopula,
               simulate, tune, validate, fit
using DataFrames
using CSV
using Dates
using LinearAlgebra
using Statistics
using StatsBase
using Distributions
using Distances
using HypothesisTests
using Plots
using StatsPlots
using Colors
using JLD2
using FileIO
using PrettyTables
using Printf
using Random
using TOML

Random.seed!(1234)

include(joinpath(_PATH_TO_SRC, "Composers.jl"))
include(joinpath(_PATH_TO_SRC, "Metrics.jl"))
include(joinpath(_PATH_TO_SRC, "Pipeline.jl"))
include(joinpath(_PATH_TO_SRC, "SyntheticMarket.jl"))
include(joinpath(_PATH_TO_SRC, "VaRBacktest.jl"))
