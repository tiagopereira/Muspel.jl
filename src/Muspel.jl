module Muspel

export Atmosphere, AtomicLine, AtomicContinuum, AtomicModel
export read_atom
export saha_boltzmann
export α_atoms_bf, σ_atoms_bf, σ_atoms_bf_tables
export α_cont, Tables_σ

using AtomicData
using Interpolations
using PeriodicTable
using StaticArrays
using Transparency
using Unitful
using YAML
import PhysicalConstants.CODATA2018: h, k_B, R_∞, c_0, m_e, m_u, e, ε_0, a_0
using ProgressMeter

@derived_dimension NumberDensity Unitful.𝐋^-3

include("types.jl")
include("read_utils.jl")
include("lte.jl")
include("background.jl")

end
