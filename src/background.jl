"""
Functions for computing background extinction.
"""

const ABUNDANCES = get_solar_abundances()
const HC_K = ustrip((h * c_0 / k_B) |> u"K * nm")


#=----------------------------------------------------------------------------
                        Cross-section tables.
----------------------------------------------------------------------------=#


"""
    function Tables_σ(
        λ::AbstractVector{T},
        log_temp::AbstractRange{T},
        log_ne::AbstractRange{T},
        atoms::AbstractVector{AtomicModel}
    ) where T <: AbstractFloat

Structure for the cross-section interpolation tables used for the background extinction.
Each field is an array of CubicSplineInterpolation functions. There is one function per
wavelength.

# Arguments:
- `λ`: Wavelengths for which tables are computed.
- `log_temp`: Range of log10 temperatures for the tables.
- `log_ne`: Range of log10 electron densities for the tables.
- `atoms`: Background atoms with continua to include.

# Fields:
- `table_nh::Vector{Interpolations.Extrapolation}`:
            Takes arguments (log10 temperature, log10 electron density). Contains
            contribution from atomic continua, H- ff and H- bf. Multiply by hydrogen
            density to get extinction per m^-1.
- `table_ne::Vector{Interpolations.Extrapolation}`:
            Takes argument (log10 temperature). Contains hydrogenic ff contribution from
            protons. Multiply by proton and electron density to get extinction per m^-1.
"""
struct Tables_σ{T0 <: Vector{Interpolations.Extrapolation}}
    table_nh::T0
    table_ne::T0
    function Tables_σ(
        λ::AbstractVector{T},
        log_temp::AbstractRange{T},
        log_ne::AbstractRange{T},
        atoms::AbstractVector{AtomicModel}
    ) where T <: AbstractFloat
        σ_atom_tables = σ_atoms_bf_tables(atoms)
        λ_u = λ * u"nm"
        temperature = 10 .^ log_temp
        temperature_u = collect(temperature) * u"K"
        electron_density = 10 .^ log_ne
        nλ = length(λ)
        nT = length(log_temp)
        ne = length(log_ne)
        T1 = eltype(λ)
        table_nh = Array{T1}(undef, nλ, nT, ne)
        table_ne = Array{T1}(undef, nλ, nT)
        p = ProgressMeter.Progress(ne * nT, desc="Computing table_nh ")
        Threads.@threads for index in CartesianIndices((ne, nT))
            ie, iT = Tuple(index)
            @. table_nh[:, iT, ie] = ustrip(σ_hminus_ff(λ_u, temperature_u[iT]) |> u"m^5")
            @. table_nh[:, iT, ie] += ustrip(σ_hminus_bf(λ_u, temperature_u[iT]) |> u"m^5")
            @. table_nh[:, iT, ie] *= electron_density[ie]
            table_nh[:, iT, ie] .+= σ_atoms_bf.(
                        Ref(σ_atom_tables),
                        Ref(atoms),
                        λ,
                        temperature[iT],
                        electron_density[ie]
            )
            ProgressMeter.next!(p)
        end
        ν_u = c_0 ./ λ_u
        for iT in 1:nT
            @. table_ne[:, iT] = ustrip(
                σ_hydrogenic_ff(ν_u, temperature_u[iT], 1) |> u"m^5"
            )
        end
        table_nh = [CubicSplineInterpolation(
                          (log_temp, log_ne), table_nh[iλ, :, :], extrapolation_bc=Line()
                   ) for iλ in 1:nλ]
        table_ne = [CubicSplineInterpolation(
                          (log_temp), table_ne[iλ, :], extrapolation_bc=Line()
                   ) for iλ in 1:nλ]
        return new{Vector{Interpolations.Extrapolation}}(table_nh, table_ne)
    end
end



#=----------------------------------------------------------------------------
            Functions for computing the background extinction.
----------------------------------------------------------------------------=#

"""
    function α_cont_no_atoms(
        λ::T,
        temperature::T,
        electron_density::T,
        h_ground_density::T,
        h_neutral_density::T,
        proton_density::T
    ) where T <: AbstractFloat

Calculates continuum extinction according to recipe in RH using functions from Transparecy.

# Arguments
- `λ`: Wavelength in nm.
- `temperature`: Temperature in K.
- `electron_density`: Electron number density in m^-3.
- `h_ground_density`: Number density of hydrogen in the ground state, unit m^-3.
- `h_neutral_density`: Number density of neutral hydrogen, unit m^-3.
- `proton_density`: Proton number density in m^-3.

# Returns
- `α`: Continuous extinction (Float) in m^-1.
"""
function α_cont_no_atoms(
    λ::T,
    temperature::T,
    electron_density::T,
    h_ground_density::T,
    h_neutral_density::T,
    proton_density::T
) where T <: AbstractFloat
    λ *= u"nm"
    temperature *= u"K"
    electron_density *= u"m^-3"
    h_ground_density *= u"m^-3"
    h_neutral_density *= u"m^-3"
    proton_density *= u"m^-3"
    α = α_hminus_ff(λ, temperature, h_neutral_density,  electron_density)
    α += α_hminus_bf(λ, temperature, h_neutral_density, electron_density)
    α += α_hydrogenic_ff(c_0 / λ, temperature, electron_density, proton_density, 1)
    α += α_h2plus_ff(λ, temperature, h_neutral_density, proton_density)
    α += α_h2plus_bf(λ, temperature, h_neutral_density, proton_density)
    α += α_thomson(electron_density)
    α += α_rayleigh_h(λ, h_ground_density)
    return ustrip(α |> u"m^-1")
end

"""
function α_cont(
    atoms::AbstractVector{AtomicModel},
    σ_atom_tables::Vector{Vector{Interpolations.FilledExtrapolation}},
    λ::T,
    temperature::T,
    electron_density::T,
    h_ground_density::T,
    h_neutral_density::T,
    proton_density::T
) where T <: AbstractFloat


Calculates continuum extinction according to recipe in RH using functions from Transparecy.
Includes continua from atom files. Atoms are treated in LTE.

# Arguments
- `atoms`: Atoms with continua to include.
- `σ_atom_tables`: Interpolation functions from σ_atoms_bf_tables().
- `λ`: Wavelength in nm.
- `temperature`: Temperature in K.
- `electron_density`: Electron number density in m^-3.
- `h_ground_density`: Number density of hydrogen in the ground state, unit m^-3.
- `h_neutral_density`: Number density of neutral hydrogen, unit m^-3.
- `proton_density`: Proton number density in m^-3.

# Returns
- `α`: Continuous extinction (Float) in m^-1.
"""
function α_cont(
    atoms::AbstractVector{AtomicModel},
    σ_atom_tables::Vector{Vector{Interpolations.FilledExtrapolation}},
    λ::T,
    temperature::T,
    electron_density::T,
    h_ground_density::T,
    h_neutral_density::T,
    proton_density::T
) where T <: AbstractFloat
    α = α_atoms_bf(
        σ_atom_tables, atoms, λ, temperature, electron_density, h_neutral_density)
    α += α_cont_no_atoms(
        λ,
        temperature,
        electron_density,
        h_ground_density,
        h_neutral_density,
        proton_density
    )
    return α
end

"""
    function α_cont_fromtables(
        tables::Tables_σ,
        iλ::Integer,
        λ::T,
        temperature::T,
        electron_density::T,
        h_ground_density::T,
        h_neutral_density::T,
        proton_density::T
    ) where T<: AbstractFloat

Computes the extinction using interpolation tables of the cross-sections. This is
significally faster than the non-tabulated version when the computation includes atomic
continua.

# Arguments
- `tables`: Pre-computed interpolation functions for "cross-sections".
- `iλ`: The index in `tables` for the wavelength λ.
- `λ`: Wavelength in nm.
- `temperature`: Temperature in K.
- `electron_density`: Electron number density in m^-3.
- `h_ground_density`: Number density of hydrogen in the ground state, unit m^-3.
- `h_neutral_density`: Number density of neutral hydrogen, unit m^-3.
- `proton_density`: Proton number density in m^-3.

# Returns
- `α`: Continuous extinction (Float) in m^-1.
"""
function α_cont_fromtables(
    tables::Tables_σ,
    iλ::Integer,
    λ::T,
    temperature::T,
    electron_density::T,
    h_ground_density::T,
    h_neutral_density::T,
    proton_density::T
) where T<: AbstractFloat
    log_temp = log10(temperature)
    log_ne = log10(electron_density)
    α = tables.table_nh[iλ](log_temp, log_ne) * h_neutral_density
    α += tables.table_ne[iλ](log_temp) * electron_density * proton_density
    λ *= u"nm"
    temperature *= u"K"
    electron_density *= u"m^-3"
    h_ground_density *= u"m^-3"
    h_neutral_density *= u"m^-3"
    proton_density *= u"m^-3"
    α_u = α_h2plus_ff(λ, temperature, h_neutral_density, proton_density)
    α_u += α_h2plus_bf(λ, temperature, h_neutral_density, proton_density)
    α_u += α_thomson(electron_density)
    α_u += α_rayleigh_h(λ, h_ground_density)
    α += ustrip(α_u |> u"m^-1")
    return α
end


#=----------------------------------------------------------------------------
            Functions to compute σ/α from background atoms.
----------------------------------------------------------------------------=#

"""
    function σ_atoms_bf_tables(atoms::AbstractVector{AtomicModel})

Returns interpolation functions for the bound-free cross section data multiplied with
abundances for each atom.

# Arguments
- `atoms`: A Vector of AtomicModels with continua.

# Returns
- `σ_atom_tables::Vector{Vector{Interpolations.FilledExtrapolation}}`:
    Interpolation functions.
"""
function σ_atoms_bf_tables(atoms::AbstractVector{AtomicModel})
    σ_atom_tables = [
        Vector{Interpolations.FilledExtrapolation}(
            undef, length(atom.continua)) for atom in atoms
    ]
    for (i, atom) in enumerate(atoms)
        for j in 1:length(atom.continua)
            σ_atom_tables[i][j] = LinearInterpolation(
                atom.continua[j].λ,
                atom.continua[j].σ * ABUNDANCES[atom.element],
                extrapolation_bc=0
            )
        end
    end
    return σ_atom_tables
end


"""
    function σ_atoms_bf(
        σ_atom_tables::Vector{Vector{Interpolations.FilledExtrapolation}},
        atoms::AbstractVector{AtomicModel},
        λ::T,
        temperature::T,
        electron_density::T
    ) where T <: AbstractFloat

Computes the bound-free cross-sections from atom files multiplied with the population of
the lower level of each bound-free transition relative to the total hydrogen population.
To get total bound-free extinction multiply with hydrogen_density.

# Arguments
- `σ_atom_tables`: Interpolation functions from σ_atoms_bf_tables().
- `atoms`: A vector of AtomicModels with continua.
- `λ`: Wavelength in nm.
- `temperature`: Temperature in K.
- `electron_density`: Number density in m^-3.

# Returns
- `σ_λ`: Sum of all cross-sections multiplied with relative populations, abundancies and
    correction for stimulated emission (Float).
"""
function σ_atoms_bf(
    σ_atom_tables::Vector{Vector{Interpolations.FilledExtrapolation}},
    atoms::AbstractVector{AtomicModel},
    λ::T,
    temperature::T,
    electron_density::T
) where T <: AbstractFloat
    σ_λ = 0
    for (i, atom) in enumerate(atoms)
        populations = saha_boltzmann(atom, temperature, electron_density)
        for j in 1:length(atom.continua)
            σ_j = σ_atom_tables[i][j](λ)
            σ_λ += (σ_j * populations[atom.continua[j].lo])
        end
    end
    return σ_λ * (1 - exp(-HC_K / (λ * temperature)))
end

"""
    function α_atoms_bf(
        σ_atom_tables::Vector{Vector{Interpolations.FilledExtrapolation}},
        atoms::AbstractVector{AtomicModel},
        λ::T,
        temperature::T,
        electron_density::T,
        hydrogen_density::T
    ) where T <: AbstractFloat


Computes the extinction per meter from continua in atoms.

# Arguments
- `σ_atom_tables`: Interpolation tables from σ_atoms_bf_tables().
- `atoms`: Vector of AtomicModels with continua.
- `λ` : Wavelength in nm.
- `temperature`: temperature in K.
- `electron_density`: Electron density in m^-3.
- `hydrogen_density`: Total hydrogen density in m^-3.

# Returns
- `α_λ`: Extinction per meter from bound-free transitions (Float).
"""
function α_atoms_bf(
    σ_atom_tables::Vector{Vector{Interpolations.FilledExtrapolation}},
    atoms::AbstractVector{AtomicModel},
    λ::T,
    temperature::T,
    electron_density::T,
    hydrogen_density::T
) where T <: AbstractFloat
    σ_λ = σ_atoms_bf(σ_atom_tables, atoms, λ, temperature, electron_density)
    return σ_λ * hydrogen_density
end
