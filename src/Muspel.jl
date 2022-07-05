module Muspel

export Atmosphere, AtomicLine, AtomicContinuum, AtomicModel
export read_atom
export saha_boltzmann

using PeriodicTable
using StaticArrays
using Transparency
using Unitful
using YAML
import PhysicalConstants.CODATA2018: h, k_B, R_∞, c_0, m_e, m_u, e, ε_0, a_0


include("types.jl")
include("read_utils.jl")
include("lte.jl")

end
