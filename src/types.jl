"""
Collection of types.
"""

abstract type AbstractAtmos{T <: Real} end

abstract type AbstractAtmos1D{N, T <: Real} <: AbstractAtmos{T} end

abstract type AbstractAtmos3D{T <: Real} <: AbstractAtmos{T} end


"""
Type for 1D atmospheres. Can contain both 1D only, or 1.5D atmospheres (multiple
columns of 1D atmospheres).
"""
struct Atmosphere1D{
    N,
    T <: Real,
    A <: AbstractArray{T, N},
    V <:AbstractVector{T},
} <: AbstractAtmos1D{N, T}
    nx::Int64
    ny::Int64
    nz::Int64
    z::V
    temperature::A
    velocity_z::A
    electron_density::A
    hydrogen1_density::A  # neutral hydrogen across all levels
    proton_density::A
end


struct Atmosphere3D{
    T <: Real,
    A <: AbstractArray{T, 3},
    V <:AbstractVector{T}
} <: AbstractAtmos3D{T}
    nx::Int64
    ny::Int64
    nz::Int64
    x::V
    y::V
    z::V
    temperature::A
    velocity_x::A
    velocity_y::A
    velocity_z::A
    electron_density::A
    hydrogen1_density::A  # neutral hydrogen across all levels
    proton_density::A
end


Base.ndims(a::AbstractAtmos) = ndims(a.temperature)


Base.size(a::AbstractAtmos) = size(a.temperature)


function Base.getindex(a::AbstractAtmos1D, args...)
    nD = length(args)
    if nD != ndims(a)
        throw(ArgumentError("Invalid number of arguments. Expected $(ndims(a)), got $nD."))
    end
    indices = to_indices(a.temperature, args)
    nx = 1
    ny = 1
    if nD == 1
        nz = length(indices[1])
    elseif nD == 2
        nz = length(indices[1])
        ny = length(indices[2])
    elseif nD == 3
        nz = length(indices[1])
        ny = length(indices[2])
        nx = length(indices[3])
    end
    if (nx == 0) | (ny == 0) | (nz == 0)
        throw(ArgumentError("All slices must have non-zero length"))
    end
    if nz == 1  # horizontal slice unsupported because velocities are not available
        throw(ArgumentError("Unsupported slice of Atmosphere1D with single z value"))
    end
    return Atmosphere1D(
        nx,
        ny,
        nz,
        a.z[indices[1]],
        a.temperature[indices...],
        a.velocity_z[indices...],
        a.electron_density[indices...],
        a.hydrogen1_density[indices...],
        a.proton_density[indices...],
    )
end


function Base.getindex(a::AbstractAtmos3D, i, j, k)
    iz, iy, ix = to_indices(a.temperature, (i, j, k))
    nx = length(ix)
    ny = length(iy)
    nz = length(iz)
    if (nx == 0) | (ny == 0) | (nz == 0)
        throw(ArgumentError("All slices must have non-zero length"))
    end
    if (nz > 1) & (ny > 1) & (nx > 1)  # simple slice, return same type
        return Atmosphere3D(
            nx,
            ny,
            nz,
            a.x[ix],
            a.y[iy],
            a.z[iz],
            a.temperature[iz, iy, ix],
            a.velocity_x[iz, iy, ix],
            a.velocity_y[iz, iy, ix],
            a.velocity_z[iz, iy, ix],
            a.electron_density[iz, iy, ix],
            a.hydrogen1_density[iz, iy, ix],
            a.proton_density[iz, iy, ix],
        )
    else  # reduced dimensionality slice, return Atmosphere1D
        if length(iz) == 1  # special case of horizontal slice
            if length(iy) > 1  # first case, swap z for y axis
                nx = length(ix)
                ny = length(iz)
                nz = length(iy)
                z = a.y[iy]
                vz = a.velocity_y[iz, iy, ix]
            elseif length(iy) == 1  # second case, swap z for x axis
                nx = length(iz)
                ny = length(iy)
                nz = length(ix)
                z = a.x[ix]
                vz = a.velocity_x[iz, iy, ix]
            end
        else
            nx = length(ix)
            ny = length(iy)
            nz = length(iz)
            z = a.z[iz]
            vz = a.velocity_z[iz, iy, ix]
        end
        return Atmosphere1D(
            nx,
            ny,
            nz,
            z,
            a.temperature[iz, iy, ix],
            vz,
            a.electron_density[iz, iy, ix],
            a.hydrogen1_density[iz, iy, ix],
            a.proton_density[iz, iy, ix],
        )
    end
end


abstract type AbstractBroadening{T <: AbstractFloat} end


struct LineBroadening{N, T} <: AbstractBroadening{T}
    natural::T
    coeff::SVector{N, T}
    temp_exp::SVector{N, T}
    hydrogen_exp::SVector{N, T}
    electron_exp::SVector{N, T}
end


struct AtomicLine{
    N,
    FloatT <: AbstractFloat,
    IntT <: Integer,
    Vl <:AbstractVector{FloatT},
    Vb <:AbstractVector{<: Real},
    Vp <:AbstractVector{<: Real},
    Vr <:AbstractVector{<: Real},
    S <: Union{Nothing, String},
}
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
    λ::Vl
    PRD::Bool
    Voigt::Bool
    label_up::S
    label_lo::S
    γ::LineBroadening{N, FloatT}
    σr_strength::Vr
    σr_shift::Vr
    π_strength::Vp
    π_shift::Vp
    σb_strength::Vb
    σb_shift::Vb
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
    τ_one_height::Vector{T}
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
        τ_one_height = Vector{T}(undef, nλ)
        source_function = Vector{T}(undef, ndep)
        α_total = Vector{T}(undef, ndep)
        α_c = Vector{T}(undef, ndep)
        j_c = Vector{T}(undef, ndep)
        ΔλD = Vector{T}(undef, ndep)
        γ = Vector{T}(undef, ndep)
        int_tmp = Vector{T}(undef, ndep)
        new{T}(
            ndep,
            nλ,
            intensity,
            τ_one_height,
            source_function,
            α_total,
            α_c,
            j_c,
            ΔλD,
            γ,
            int_tmp
        )
    end
end


struct RTBufferLTE{T <: AbstractFloat}
    ndep::Int
    nλ::Int
    intensity::Vector{T}
    source_function::Vector{T}
    α_total::Vector{T}
    α_l::Vector{T}
    α_c::Vector{T}
    ΔλD::Vector{T}
    γ::Vector{T}
    int_tmp::Vector{T}
    function RTBufferLTE(ndep, nλ, f_type)
        T = f_type
        intensity = Vector{T}(undef, nλ)
        source_function = Vector{T}(undef, ndep)
        α_total = Vector{T}(undef, ndep)
        α_l = Vector{T}(undef, ndep)
        α_c = Vector{T}(undef, ndep)
        ΔλD = Vector{T}(undef, ndep)
        γ = Vector{T}(undef, ndep)
        int_tmp = Vector{T}(undef, ndep)
        new{T}(
            ndep,
            nλ,
            intensity,
            source_function,
            α_total,
            α_l,
            α_c,
            ΔλD,
            γ,
            int_tmp
        )
    end
end


struct RTBufferStokes{T <: AbstractFloat}
    ndep::Int
    nλ::Int
    stokes::Array{T, 2}
    profiles::Array{T, 2}
    α_c::Vector{T}
    j_c::Vector{T}
    α_l::Vector{T}
    j_l::Vector{T}
    αI::Vector{T}
    ΔλD::Vector{T}
    adamp::Vector{T}
    int_tmp::Array{T, 2}
    function RTBufferStokes(ndep, nλ; t::Type{T}=Float32) where T
        stokes = Array{T}(undef, nλ, 4)
        profiles = Array{T}(undef, 6, ndep)
        α_c = Vector{T}(undef, ndep)
        j_c = Vector{T}(undef, ndep)
        α_l = Vector{T}(undef, ndep)
        j_l = Vector{T}(undef, ndep)
        αI = Vector{T}(undef, ndep)
        ΔλD = Vector{T}(undef, ndep)
        adamp = Vector{T}(undef, ndep)
        int_tmp = Array{T}(undef, ndep, 4)
        new{T}(ndep, nλ, stokes, profiles, α_c, j_c, α_l, j_l, αI, ΔλD, adamp, int_tmp)
    end
end
