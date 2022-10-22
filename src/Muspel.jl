module Muspel

export Atmosphere, AtomicLine, AtomicContinuum, AtomicModel
export AbstractBroadening, LineBroadening
export read_atom
export saha_boltzmann, saha_boltzmann!
export AbstractExtinctionItp, ExtinctionItpLTE, ExtinctionItpNLTE
export α_cont, α_cont_no_itp
export get_atoms_bf_interpolant, create_σ_itp_LTE, create_σ_itp_NLTE
export piecewise_1D_nn, piecewise_1D_linear

using AtomicData
using Interpolations
using PeriodicTable
using StaticArrays
using Transparency
using Unitful
using YAML
import PhysicalConstants.CODATA2018: h, k_B, R_∞, c_0, m_e, m_u, e, ε_0, a_0
using ProgressMeter

include("types.jl")
include("read_utils.jl")
include("lte.jl")
include("background.jl")

end
