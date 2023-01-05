"""
Calculate quantities in local thermodynamical equilibrium (LTE).
"""

const saha_const_u = ustrip(h^2 / (2 * π * m_e * k_B))
const k_B_u = ustrip(k_B)
const Hχ∞ = 2.1787094174620437e-18  # in J, from NIST retrieved Jan 2023


"""
    function saha_boltzmann(
        χ::SVector,
        g::SVector,
        stage::SVector,
        temperature::T,
        electron_density::T,
        atom_density::T
    ) where T <: AbstractFloat

    function saha_boltzmann(
        atom::AtomicModel,
        temperature::T,
        electron_density::T,
        atom_density::T,
    ) where T <: AbstractFloat

Calculates atomic level populations according to the Saha-Boltzmann
distribution.

# Arguments
- `χ`: level energies in J
- `g`: statistical weights of levels
- `stage`: ionisation stage of each level (starting at 1 for neutral)
- Or, instead of the three above, an instance of `AtomicModel`
- `temperature`: temperature in Kelvin
- `electron_density`: electron density in m^-3
- `atom_density`: total number density (in all levels) of target species

# Returns
- `populations`: MVector{nlevels} with relative level populations in m^-3
"""
function saha_boltzmann(
    χ::SVector{N, <: Real},
    g::SVector{N, <: Real},
    stage::SVector{N, <: Real},
    temperature::T,
    electron_density::T,
    atom_density::Real,
) where {N, T <: AbstractFloat}
    populations = MVector{length(χ), T}(undef)
    saha_boltzmann!(χ, g, stage, temperature, electron_density, atom_density, populations)
    return populations
end

function saha_boltzmann(
    atom::AtomicModel,
    temperature::T,
    electron_density::T,
    atom_density::Real,
) where T <: AbstractFloat
    saha_boltzmann(atom.χ, atom.g, atom.stage, temperature, electron_density, atom_density)
end


"""
    function saha_boltzmann!(
        χ::SVector,
        g::SVector,
        stage::SVector,
        temperature::T,
        electron_density::T,
        atom_density::T,
        populations::AbstractArray{T, 1},
    ) where T <: AbstractFloat

    function saha_boltzmann!(
        atom::AtomicModel,
        temperature::T,
        electron_density::T,
        atom_density::T,
        populations::AbstractArray{T, 1},
    ) where T <: AbstractFloat

Inplace version of `saha_boltzmann``. Calculates atomic level populations
according to the Saha-Boltzmann distribution, placing them in an existing
`populations` array.

# Arguments
- `χ`: level energies in J
- `g`: statistical weights of levels
- `stage`: ionisation stage of each level (starting at 1 for neutral)
- Or, instead of the three above, an instance of `AtomicModel`
- `temperature`: temperature in Kelvin
- `electron_density`: electron density in m^-3
- `atom_density`: total number density (in all levels) of target species
- `populations`: 1D array for output, must be same length as number of levels
"""
function saha_boltzmann!(
    χ::SVector{N, <: Real},
    g::SVector{N, <: Real},
    stage::SVector{N, <: Real},
    temperature::Real,
    electron_density::Real,
    atom_density::Real,
    populations::AbstractVector{T},
) where {N, T <: AbstractFloat}
    nlevels = length(χ)
    @assert nlevels == length(populations)
    kT = convert(T, k_B_u * temperature)
    saha_factor = convert(T, (saha_const_u / temperature) ^ (3/2) * electron_density / 2)
    total = one(T)
    for i in 2:nlevels
        ΔE = χ[i] - χ[1]
        ΔZ = stage[i] - stage[1]
        populations[i] = g[i] / g[1] * exp(-ΔE / kT)
        for s in 1:ΔZ
            populations[i] /= saha_factor
        end
        total += populations[i]
    end
    populations[1] = 1
    for i = 1:nlevels
        populations[i] *= atom_density / total
    end
    return nothing
end


function saha_boltzmann!(
    atom::AtomicModel,
    temperature::Real,
    electron_density::Real,
    atom_density::Real,
    populations::AbstractVector{<: Real},
)
    saha_boltzmann!(
        atom.χ,
        atom.g,
        atom.stage,
        temperature,
        electron_density,
        atom_density,
        populations
    )
    return nothing
end


"""
    h_saha(temp::T, electron_density::T)::T where {T <: Real}

Calculate ionisation fraction of hydrogen using Saha.
"""
function h_ionfrac_saha(temp::T, electron_density::T)::T where {T <: Real}
    saha = (temp / saha_const_u) ^ (3/2) / electron_density * exp(-Hχ∞ / (k_B_u * temp))
    return 1 - (1 / (1 + saha))
end
