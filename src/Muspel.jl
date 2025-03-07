module Muspel

export Atmosphere1D, Atmosphere3D
export AtomicLine, AtomicContinuum, AtomicModel
export RTBuffer
export AbstractBroadening, LineBroadening
export read_atom
export saha_boltzmann, saha_boltzmann!, create_nstar_itp
export AbstractExtinctionItp, ExtinctionItpLTE, ExtinctionItpNLTE
export α_cont, α_cont_no_itp, get_σ_itp
export piecewise_1D_nn, piecewise_1D_linear, piecewise_1D_linear!, piecewise_1D_bezier3!,
       feautrier
export read_atmos_rh, read_atmos_rh_index, read_atmos_hpops_rh, read_pops_rh
export read_atmos_hpops_multi3d, read_atmos_multi3d, read_pops_multi3d
export doppler_width, damping, calc_broadening, create_voigt_itp
export blackbody_λ
export calc_line_1D!, calc_τ_cont!, calc_τ_one_height
export incline_atmos, incline_data!, incline_data_inv!, project_vector!
export g_eff, g_LS, g_jj, g_JK, zeeman_strength, get_zeeman_components

using AtomicData
using Interpolations
using PeriodicTable
using StaticArrays
using Transparency
using Unitful
using YAML
import PhysicalConstants.CODATA2018: h, k_B, R_∞, c_0, m_e, m_u, e, ε_0, a_0
import SpecialFunctions: erfcx
using ProgressMeter


include("types.jl")
include("read_utils.jl")
include("lte.jl")
include("line.jl")
include("background.jl")
include("formal_solvers.jl")
include("read_atmos.jl")
include("utils.jl")
include("intensity.jl")
include("incline.jl")
include("zeeman.jl")

end
