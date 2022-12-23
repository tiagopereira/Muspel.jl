"""
Collection of types.
"""

abstract type AbstractAtmos{T <: AbstractFloat} end

struct Atmosphere{FloatT <: AbstractFloat}
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
    function Atmosphere(
        x::AbstractArray{FloatT, 1},
        y::AbstractArray{FloatT, 1},
        z::AbstractArray{FloatT, 1},
        temperature::AbstractArray{FloatT, 3},
        velocity_z::AbstractArray{FloatT, 3},
        electron_density::AbstractArray{FloatT, 3},
        hydrogen_density::AbstractArray{FloatT, 4},
    ) where FloatT <: AbstractFloat
        nz, ny, nx, nh_levels = size(hydrogen_density)
        @assert (nz, ny, nx) == (length(z), length(y), length(x))
        @assert size(temperature) == (nz, ny, nx)
        @assert size(velocity_z) == (nz, ny, nx)
        @assert size(electron_density) == (nz, ny, nx)
        if nh_levels == 1  # hydrogen populations not given, must use Saha to get ionisation
            hpops = view(hydrogen_density, :, :, :, 1)
            ATOM_FILE = "H_3.yaml"
            filepath = joinpath(@__DIR__, "..", "data", "atoms", ATOM_FILE)
            if isfile(filepath)
                H_atom = read_atom(filepath)
                hydrogen1_density = Array{FloatT}(undef, nz, ny, nx)
                proton_density = similar(hydrogen1_density)
                tmp = Vector{FloatT}(undef, H_atom.nlevels)

                for i in eachindex(temperature)
                    saha_boltzmann!(
                        H_atom,
                        temperature[i],
                        electron_density[i],
                        hpops[i],
                        tmp
                    )
                    proton_density[i] = tmp[end]
                    hydrogen1_density[i] = max(0, hpops[i] - proton_density[i])
                end
            else
                error("Hydrogen model atom $ATOM_FILE was not found.")
            end
        else
            hydrogen1_density = sum(hydrogen_density[:, :, :, 1:end-1], dims=4)[:, :, :, 1]
            proton_density = hydrogen_density[:, :, :, end]
        end
        new{FloatT}(nx, ny, nz,
                    x, y, z,
                    temperature, velocity_z, electron_density,
                    hydrogen1_density, proton_density)
    end
end


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
