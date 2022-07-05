"""
Calculate quantities in local thermodynamical equilibrium (LTE).
"""

const saha_const_u = ustrip(h^2 / (2 * π * m_e * k_B))
const k_B_u = ustrip(k_B)


"""
    function saha_boltzmann(χ::SVector,
                            g::SVector,
                            stage::SVector,
                            temperature::T,
                            electron_density::T) where T <: AbstractFloat

    function saha_boltzmann(atom::AtomicModel,
                            temperature::T,
                            electron_density::T) where T <: AbstractFloat

Calculates atomic level populations according to the Saha-Boltzmann
distribution.

# Arguments
- `χ`: level energies in J
- `g`: statistical weights of levels
- `stage`: ionisation stage of each level (starting at 1 for neutral)
- Or, instead of the three above, an instance of `AtomicModel`
- `temperature`: temperature in Kelvin
- `electron_density`: electron density in m^-3

# Returns
- `populations`: MVector{nlevels} with relative level populations in m^-3
"""
function saha_boltzmann(χ::SVector,
                        g::SVector,
                        stage::SVector,
                        temperature::T,
                        electron_density::T) where T <: AbstractFloat
    nlevels = length(χ)
    populations = MVector{nlevels, T}(undef)
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
    populations / total
end

function saha_boltzmann(atom::AtomicModel,
                        temperature::T,
                        electron_density::T) where T <: AbstractFloat
    saha_boltzmann(atom.χ, atom.g, atom.stage, temperature, electron_density)
end
