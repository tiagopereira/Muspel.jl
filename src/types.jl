"""
Collection of types.
"""

abstract type AbstractAtmos{T <: AbstractFloat} end

struct Atmosphere{FloatT <: AbstractFloat} <: AbstractAtmos{FloatT}
    nx::Int64
    ny::Int64
    nz::Int64
    x::Array{FloatT, 1}
    y::Array{FloatT, 1}
    z::Array{FloatT, 1}
    temperature::Array{FloatT, 3}
    velocity_z::Array{FloatT, 3}
    electron_density::Array{FloatT, 3}
    hydrogen1_density::Array{FloatT, 3}  # neutral hydrogen across all levels
    proton_density::Array{FloatT, 3}
end


# Basic Multi3D atmos using native shape. Does not calculate Saha-Boltzmann.
struct AtmosphereM3D{FloatT <: AbstractFloat} <: AbstractAtmos{FloatT}
    nx::Int64
    ny::Int64
    nz::Int64
    z::Array{FloatT, 1}
    temperature::Array{FloatT, 3}
    velocity_z::Array{FloatT, 3}
    electron_density::Array{FloatT, 3}
    hydrogen1_density::Array{FloatT, 3}  # neutral hydrogen across all levels
    proton_density::Array{FloatT, 3}
end


struct Atmosphere1D{FloatT <: AbstractFloat} <: AbstractAtmos{FloatT}
    nz::Int64
    z::Array{FloatT, 1}
    temperature::Array{FloatT, 1}
    velocity_z::Array{FloatT, 1}
    electron_density::Array{FloatT, 1}
    hydrogen1_density::Array{FloatT, 1}
    proton_density::Array{FloatT, 1}
end


Base.getindex(a::Atmosphere, i, j) = Atmosphere1D(
    a.nz,
    a.z,
    a.temperature[:,i, j],
    a.velocity_z[:,i, j],
    a.electron_density[:,i, j],
    a.hydrogen1_density[:,i, j],
    a.proton_density[:, i, j],
)


Base.getindex(a::AtmosphereM3D, i, j) = Atmosphere1D(
    length(a.z),
    a.z,
    a.temperature[i, j, :],
    a.velocity_z[i, j, :],
    a.electron_density[i, j, :],
    a.hydrogen1_density[i, j, :],
    a.proton_density[i, j, :],
)


abstract type AbstractBroadening{T <: AbstractFloat} end


struct LineBroadening{N, M, T} <: AbstractBroadening{T}
    natural::T
    hydrogen_const::SVector{N, T}
    hydrogen_exp::SVector{N, T}
    electron_const::SVector{M, T}
    electron_exp::SVector{M, T}
    stark_linear_const::T
    stark_linear_exp::T
end


struct AtomicLine{N, M, FloatT <: AbstractFloat, IntT <: Integer}
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
    mass::FloatT
    λ::Vector{FloatT}
    PRD::Bool
    Voigt::Bool
    label_up::String
    label_lo::String
    γ::LineBroadening{N, M, FloatT}
end



struct AtomicContinuum{Nλ, FloatT <: AbstractFloat, IntT <: Integer}
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


struct RTBuffer{T <: AbstractFloat}
    ndep::Int
    nλ::Int
    intensity::Vector{T}
    source_function::Vector{T}
    α_total::Vector{T}
    α_c::Vector{T}
    j_c::Vector{T}
    ΔλD::Vector{T}
    γ::Vector{T}
    int_tmp::Vector{T}
    function RTBuffer(ndep, nλ, f_type)
        T = f_type
        intensity = Vector{T}(undef, nλ)
        source_function = Vector{T}(undef, ndep)
        α_total = Vector{T}(undef, ndep)
        α_c = Vector{T}(undef, ndep)
        j_c = Vector{T}(undef, ndep)
        ΔλD = Vector{T}(undef, ndep)
        γ = Vector{T}(undef, ndep)
        int_tmp = Vector{T}(undef, ndep)
        new{T}(ndep, nλ, intensity, source_function, α_total, α_c, j_c, ΔλD, γ, int_tmp)
    end
end
