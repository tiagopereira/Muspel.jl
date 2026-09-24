"""
Collection of types.
"""

abstract type AbstractAtmosphere{N, T <: Real} end


"""
    Atmosphere{N, T, A, V, Vel, B}

Type for model atmospheres. The dimensionality `N` refers to the data arrays
(following the order z, y, x for the dimensions that exist). Vector quantities
are stored in NamedTuples, whose keys encode which components are present:

* `velocity`: any subset of `(:x, :y, :z)`. Atmospheres for 1D/1.5D work
  typically carry only `(z = ...,)` to save memory and load times.
* `magnetic_field`: `nothing`, or a NamedTuple with keys `(:x, :y, :z)`.

All components must have the same array type as the scalar fields
(`AbstractArray{T, N}`). Slicing preserves all components: no field is
dropped or renamed, only the axis vectors are adjusted.
"""
struct Atmosphere{
    N,
    T <: Real,
    A <: AbstractArray{T, N},
    V <: AbstractVector{T},
    Vel <: NamedTuple,
    B <: Union{Nothing, NamedTuple},
} <: AbstractAtmosphere{N, T}
    nx::Int64
    ny::Int64
    nz::Int64
    x::V
    y::V
    z::V
    temperature::A
    velocity::Vel
    magnetic_field::B
    electron_density::A
    hydrogen1_density::A  # neutral hydrogen across all levels
    proton_density::A

    function Atmosphere(
        nx::Integer,
        ny::Integer,
        nz::Integer,
        x::V,
        y::V,
        z::V,
        temperature::A,
        velocity::Vel,
        magnetic_field::B,
        electron_density::A,
        hydrogen1_density::A,
        proton_density::A,
    ) where {
        N,
        T <: Real,
        A <: AbstractArray{T, N},
        V <: AbstractVector{T},
        Vel <: NamedTuple,
        B <: Union{Nothing, NamedTuple},
    }
        _check_components(velocity, A, (:x, :y, :z), "velocity")
        _check_components(magnetic_field, A, (:x, :y, :z), "magnetic_field"; exact=true)
        new{N, T, A, V, Vel, B}(
            nx, ny, nz, x, y, z, temperature, velocity, magnetic_field,
            electron_density, hydrogen1_density, proton_density,
        )
    end
end


"""
    _check_components(nt, A, valid_keys, name; exact=false)

Ensure that all elements of `nt` are arrays of type `A`, and that its keys
are valid. If `exact=true`, the keys must match `valid_keys` exactly.
"""
function _check_components(
    nt::NamedTuple, ::Type{A}, valid_keys, name; exact=false
) where A <: AbstractArray
    ks = keys(nt)
    if exact && ks != valid_keys
        throw(ArgumentError("$name must have keys $valid_keys, got $ks"))
    elseif !all(k in valid_keys for k in ks)
        throw(ArgumentError("Invalid keys in $name: $ks. Valid keys are $valid_keys."))
    end
    for (k, v) in pairs(nt)
        if !(v isa A)
            throw(ArgumentError(
                "Component $k of $name must be of type $A, got $(typeof(v))"
            ))
        end
    end
    return nothing
end

_check_components(::Nothing, ::Type{<:AbstractArray}, valid_keys, name; exact=false) = nothing


"""
    has_magnetic_field(a::AbstractAtmosphere)

Compile-time check (resolved from the type) for the presence of a magnetic field.
"""
has_magnetic_field(a::AbstractAtmosphere) = a.magnetic_field !== nothing


Base.ndims(::AbstractAtmosphere{N}) where N = N


Base.size(a::AbstractAtmosphere) = size(a.temperature)


@inline _slice_axis(ax::AbstractVector, i::Integer) = ax[i:i]
@inline _slice_axis(ax::AbstractVector, i) = ax[i]

@inline _slice_components(nt::NamedTuple, idx) = map(v -> v[idx...], nt)
_slice_components(::Nothing, idx) = nothing


function Base.getindex(a::AbstractAtmosphere{N}, args...) where N
    nD = length(args)
    if nD != N
        throw(ArgumentError("Invalid number of arguments. Expected $N, got $nD."))
    end
    indices = to_indices(a.temperature, args)
    nz = length(indices[1])
    ny = N > 1 ? length(indices[2]) : 1
    nx = N > 2 ? length(indices[3]) : 1
    if (nx == 0) | (ny == 0) | (nz == 0)
        throw(ArgumentError("All slices must have non-zero length"))
    end
    x = a.x
    y = a.y
    z = a.z
    if N == 3
        if !isempty(a.x)
            x = _slice_axis(a.x, indices[3])
        end
        if !isempty(a.y)
            y = _slice_axis(a.y, indices[2])
        end
        z = _slice_axis(a.z, indices[1])
        if nz == 1  # horizontal slice: rotate the axis vectors, data arrays untouched
            if isempty(a.x) | isempty(a.y)
                throw(ArgumentError("Unsupported slice with a single z value"))
            elseif ny > 1  # new vertical is the y axis
                z = _slice_axis(a.y, indices[2])
                y = _slice_axis(a.z, indices[1])
                nz, ny = ny, nz
            elseif nx > 1  # new vertical is the x axis
                z = _slice_axis(a.x, indices[3])
                x = _slice_axis(a.z, indices[1])
                nz, nx = nx, nz
            else
                throw(ArgumentError("Cannot slice a single point"))
            end
        end
    else
        z = _slice_axis(a.z, indices[1])
        if N > 1 && !isempty(a.y)
            y = _slice_axis(a.y, indices[2])
        end
        if nz == 1
            throw(ArgumentError("Unsupported slice with a single z value"))
        end
    end
    return Atmosphere(
        nx,
        ny,
        nz,
        x,
        y,
        z,
        a.temperature[indices...],
        _slice_components(a.velocity, indices),
        _slice_components(a.magnetic_field, indices),
        a.electron_density[indices...],
        a.hydrogen1_density[indices...],
        a.proton_density[indices...],
    )
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
