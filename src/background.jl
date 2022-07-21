"""
Functions for computing background extinction.
"""

const _ABUNDANCES = get_solar_abundances()
# (h c)/k_B in units (K nm):
const _HC_K = (ustrip(h |> u"J*s") * ustrip(c_0 |> u"nm/s") / ustrip(k_B |> u"J*K^-1"))


#=----------------------------------------------------------------------------
                        Cross-section tables.
----------------------------------------------------------------------------=#


"""
    function Tables_σ(
        λ::Vector{T}, 
        log_temp::StepRangeLen, 
        log_ne::StepRangeLen, 
        atoms::Vector{AtomicModel}
    ) where T<:AbstractFloat

Structure for the cross-section interpolation tables used for the background extinction. 
Each field is an array of CubicSplineInterpolation functions. There is one function per  
wavelength.

# Arguments:
- `λ`: Wavelengths for which tables are computed.
- `log_temp`: Range of log10 temperatures for the tables.
- `log_ne`: Range of log10 electron densities for the tables.
- `atoms`: Background atoms with continua to include.

# Fields:
- `table_a::Vector{Interpolations.Extrapolation}`: 
            Takes arguments (log10 temperature, log10 electron density). Contains 
            contribution from atomic continua, H- ff and H- bf. Multiply with hydrogen 
            density to get extinction per m^-1.
- `table_b::Vector{Interpolations.Extrapolation}`:
            Takes argument (log10 temperature). Contains hydrogenic ff contribution from 
            protons. Multiply with proton and electron density to get extinction per m^-1.
"""
struct Tables_σ{T0 <: Vector{Interpolations.Extrapolation}}
    table_a::T0
    table_b::T0 
    function Tables_σ(
        λ::Vector{T}, 
        log_temp::StepRangeLen{T}, 
        log_ne::StepRangeLen{T}, 
        atoms::Vector{AtomicModel}
    ) where T<:AbstractFloat
        σ_atom_tables = σ_atoms_bf_tables(atoms)
        λ_u = λ*u"nm"
        temperature = 10 .^ log_temp
        temperature_u = collect(temperature)*u"K"
        electron_density = 10 .^ log_ne
        nλ = length(λ)
        nT = length(log_temp)
        ne = length(log_ne)
        T1 = typeof(λ[1])
        table_a = Array{T1}(undef, nλ, nT, ne)
        table_b = Array{T1}(undef, nλ, nT)
        p = ProgressMeter.Progress(ne * nT, desc="Computing table_a ")
        Threads.@threads for ii in CartesianIndices((ne, nT))
            ie, iT = Tuple(ii)
            @. table_a[:,iT,ie] = ustrip(σ_hminus_ff(λ_u, temperature_u[iT]) |> u"m^5")
            @. table_a[:,iT,ie] += ustrip(σ_hminus_bf(λ_u, temperature_u[iT]) |> u"m^5")
            @. table_a[:,iT,ie] *= electron_density[ie]
            table_a[:,iT,ie] .+= σ_atoms_bf.(
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
            @. table_b[:,iT] = ustrip(σ_hydrogenic_ff(ν_u, temperature_u[iT], 1) |> u"m^5")
        end
        table_a = [CubicSplineInterpolation(
            (log_temp, log_ne), table_a[iλ,:,:], extrapolation_bc=Line()
            ) for iλ in 1:nλ]
        table_b = [CubicSplineInterpolation(
            (log_temp), table_b[iλ,:], extrapolation_bc=Line()
            ) for iλ in 1:nλ]
        return new{Vector{Interpolations.Extrapolation}}(table_a, table_b)
    end
end



#=----------------------------------------------------------------------------
            Functions for computing the background extinction.
----------------------------------------------------------------------------=#

"""
    function α_cont(
        λ::Unitful.Length, 
        temperature::Unitful.Temperature,
        electron_density::NumberDensity, 
        h_ground_density::NumberDensity, 
        h_neutral_density::NumberDensity, 
        proton_density::NumberDensity
    )

    function α_cont(
        λ::T, 
        temperature::T, 
        electron_density::T, 
        h_ground_density::T, 
        h_neutral_density::T, 
        proton_density::T
    ) where T <: AbstractFloat

    function α_cont(
        atoms::Vector{AtomicModel}, 
        σ_atom_tables::Vector{Vector{Interpolations.FilledExtrapolation}}, 
        λ::T, 
        temperature::T, 
        electron_density::T, 
        h_ground_density::T, 
        h_neutral_density::T, 
        proton_density::T
    ) where T <: AbstractFloat


Calculates continuum extinction according to recipe in RH using functions from Transparecy.
Optionally includes continua from atom files. Atoms are treated in LTE.

# Arguments
- `atoms`: (optional) Atoms with continua to include.
- `σ_atom_tables`: (optional) Interpolation functions from σ_atoms_bf_tables().
- `λ`: Wavelength in nm or with a Unitful.Length unit.
- `temperature`: Unit K.
- `electron_density`: Number density, unit m^-3.
- `h_ground_density`: Number density of hydrogen in the ground state.
- `h_neutral_density`: Number density of neutral hydrogen.
- `proton_density`: Number density, unit m^-3.

# Returns
- `α`: Continuous extinction (Float) in m^-1. Has units if input had units.

"""
function α_cont(
    λ::Unitful.Length, 
    temperature::Unitful.Temperature,
    electron_density::NumberDensity, 
    h_ground_density::NumberDensity, 
    h_neutral_density::NumberDensity, 
    proton_density::NumberDensity
)
    α = α_hminus_ff(λ, temperature, h_neutral_density,  electron_density)
    α += α_hminus_bf(λ, temperature, h_neutral_density, electron_density)
    α += α_hydrogenic_ff(c_0 / λ, temperature, electron_density, proton_density, 1)
    α += α_h2plus_ff(λ, temperature, h_neutral_density, proton_density)
    α += α_h2plus_bf(λ, temperature, h_neutral_density, proton_density)
    α += α_thomson(electron_density)
    α += α_rayleigh_h(λ, h_ground_density)
    return α
end

function α_cont(
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
    α = α_cont(
        λ, 
        temperature, 
        electron_density, 
        h_ground_density, 
        h_neutral_density, 
        proton_density
    )
    return ustrip(α |> u"m^-1")
end

function α_cont(
    atoms::Vector{AtomicModel}, 
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
    α += α_cont(
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
    function α_cont(
        tables::Tables_σ, 
        iλ::Integer,
        λ::T, 
        temperature::T,
        electron_density::T, 
        h_ground_density::T, 
        h_neutral_density::T, 
        proton_density::T
    ) where T<: AbstractFloat

    function α_cont(tables::Tables_σ, iλ::Integer, args...) where args <: Unitful.Quantity

If given a Tables_σ struct, the function returns the extinction computed from 
interpolation tables of the cross-sections. This is significally faster than the 
non-tabulated version when the computation includes many atomic continua.

# Arguments
    - `tables`: Pre-computed interpolation functions for "cross-sections".
    - `iλ`: The index in `tables` for the wavelength λ.
"""
function α_cont(
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
    α = tables.table_a[iλ](log_temp, log_ne) * h_neutral_density
    α += tables.table_b[iλ](log_temp) * electron_density * proton_density
    α += ustrip(α_h2plus_ff(
            λ*u"nm", 
            temperature*u"K", 
            h_neutral_density*u"m^-3", 
            proton_density*u"m^-3"
    )|>u"m^-1")
    α += ustrip(α_h2plus_bf(
        λ*u"nm", 
        temperature*u"K", 
        h_neutral_density*u"m^-3", 
        proton_density*u"m^-3"
    )|>u"m^-1")
    α += ustrip(α_thomson(electron_density*u"m^-3")|>u"m^-1")
    α += ustrip(α_rayleigh_h(λ*u"nm", h_ground_density*u"m^-3")|>u"m^-1")
    return α
end

function α_cont(
    tables::Tables_σ, 
    iλ::Integer, 
    λ::Unitful.Length, 
    temperature::Unitful.Temperature,
    electron_density::NumberDensity, 
    h_ground_density::NumberDensity, 
    h_neutral_density::NumberDensity, 
    proton_density::NumberDensity
)
    log_temp = log10(ustrip(temperature |> u"K"))
    log_ne = log10(ustrip(electron_density |> u"m^-3"))
    α = (tables.table_a[iλ](log_temp, log_ne) * u"m^2" * h_neutral_density) |> u"m^-1"
    α += (tables.table_b[iλ](log_temp)u"m^5" * electron_density * proton_density)|>u"m^-1"
    α += α_h2plus_ff(λ, temperature, h_neutral_density, proton_density)
    α += α_h2plus_bf(λ, temperature, h_neutral_density, proton_density)
    α += α_thomson(electron_density)
    α += α_rayleigh_h(λ, h_ground_density)
    return α
end


#=----------------------------------------------------------------------------
            Functions to compute σ/α from background atoms.
----------------------------------------------------------------------------=#

"""
    function σ_atoms_bf_tables(atoms::Vector{AtomicModel})

Returns interpolation functions for the bound-free cross section data multiplied with 
abundances for each atom.

# Arguments
- `atoms`: A Vector of AtomicModels with continua.

# Returns
- `σ_atom_tables::Vector{Vector{Interpolations.FilledExtrapolation}}`: 
    Interpolation functions.
"""
function σ_atoms_bf_tables(atoms::Vector{AtomicModel})
    σ_atom_tables = [
        Vector{Interpolations.FilledExtrapolation}(
            undef, length(atom.continua)) for atom in atoms
    ]
    for (i, atom) in enumerate(atoms)
        for j in 1:length(atom.continua)
            σ_atom_tables[i][j] = LinearInterpolation(
                atom.continua[j].λ,
                atom.continua[j].σ * _ABUNDANCES[atom.element],
                extrapolation_bc=0
            )
        end
    end
    return σ_atom_tables
end


"""
    function σ_atoms_bf(
        σ_atom_tables::Vector{Vector{Interpolations.FilledExtrapolation}}, 
        atoms::Vector{AtomicModel}, 
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
    atoms::Vector{AtomicModel}, 
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
    return σ_λ * (1 - exp(- _HC_K/(λ*temperature)))
end

"""
    function α_atoms_bf(
        σ_atom_tables::Vector{Vector{Interpolations.FilledExtrapolation}}, 
        atoms::Vector{AtomicModel}, 
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
    atoms::Vector{AtomicModel}, 
    λ::T, 
    temperature::T, 
    electron_density::T, 
    hydrogen_density::T
) where T <: AbstractFloat
    σ_λ = σ_atoms_bf(σ_atom_tables, atoms, λ, temperature, electron_density)
    return σ_λ * hydrogen_density
end