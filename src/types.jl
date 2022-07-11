"""
Collection of types.
"""

struct Atmosphere{FloatT <: AbstractFloat, IntT <: Integer}
    nx::IntT
    ny::IntT
    nz::IntT
    nh_levels::IntT
    x::Array{FloatT, 1}
    y::Array{FloatT, 1}
    z::Array{FloatT, 1}
    temperature::Array{FloatT, 3}
    velocity_z::Array{FloatT, 3}
    electron_density::Array{FloatT, 3}
    hydrogen_density::Array{FloatT, 4}
    function Atmosphere(x::Array{FloatT, 1},
                        y::Array{FloatT, 1},
                        z::Array{FloatT, 1},
                        temperature::Array{FloatT, 3},
                        velocity_z::Array{FloatT, 3},
                        electron_density::Array{FloatT, 3},
                        hydrogen_density::Array{FloatT, 4}) where FloatT <: AbstractFloat
        nz, ny, nx, nh_levels = size(hydrogen_density)
        IntT = typeof(nz)
        @assert size(temperature) == (nz, ny, nx)
        @assert size(velocity_z) == (nz, ny, nx)
        @assert size(electron_density) == (nz, ny, nx)
        new{FloatT, IntT}(nx, ny, nz, nh_levels,
                          x, y, z,
                          temperature, velocity_z, electron_density, hydrogen_density)
    end
end


struct AtomicLine{N, FloatT <: AbstractFloat, IntT <: Integer}
    nλ::IntT
    χup::FloatT
    χlo::FloatT
    gup::IntT
    glo::IntT
    Aul::FloatT
    Blu::FloatT
    Bul::FloatT
    λ0::FloatT  # in nm
    f_value::FloatT
    λ::Vector{FloatT}
    PRD::Bool
    Voigt::Bool
    label_up::String
    label_lo::String
    γ_vdW_const::SVector{N, FloatT}
    γ_vdW_exp::SVector{N, FloatT}
    γ_quad_stark_const::FloatT
end


struct AtomicContinuum{Nλ, FloatT <: AbstractFloat, IntT <: Integer}
    up::IntT
    lo::IntT
    nλ::IntT
    λedge::FloatT  # in nm
    σ::SVector{Nλ, FloatT}  # m^-2
    λ::SVector{Nλ, FloatT}  # nm
end


struct AtomicContinuum2{Nλ, FloatT <: AbstractFloat, IntT <: Integer}
    up::IntT
    lo::IntT
    nλ::IntT
    λedge::FloatT  # in nm
    σ::SVector{Nλ, FloatT}  # m^-2
    λ::SVector{Nλ, FloatT}  # nm
end


struct AtomicModel{Nlevel, FloatT <: AbstractFloat, IntT <: Integer}
    element::Symbol
    nlevels::IntT
    nlines::IntT
    ncontinua::IntT
    Z::IntT
    mass::FloatT
    χ::SVector{Nlevel, FloatT}  # Energy in J or aJ?
    g::SVector{Nlevel, IntT}
    stage::SVector{Nlevel, IntT}
    label::Vector{String}
    lines::Vector{AtomicLine}
    continua::Vector{AtomicContinuum}
end
