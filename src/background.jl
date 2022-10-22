"""
Functions for computing background extinction.
"""

const ABUNDANCES = get_solar_abundances()
const HC_K = ustrip((h * c_0 / k_B) |> u"K * nm")
const σ_THOMSON = ustrip(σ_thomson)


#=----------------------------------------------------------------------------
                        Types
----------------------------------------------------------------------------=#


abstract type AbstractExtinctionItp{T <: AbstractFloat} end


"""
Interpolant structure for continuum extinction, for use when hydrogen
populations are computed in LTE (simplifies calculations).
"""
struct ExtinctionItpLTE{
    T,
    ITP_2D <: Interpolations.AbstractInterpolation{T, 2}
}  <: AbstractBroadening{T}
    σ_H::ITP_2D
    σ_H2::ITP_2D
    λ::T
end


"""
Interpolant structure for continuum extinction, for use when hydrogen
populations are given explicitly (e.g. non-equilibrium ionisation or NLTE).
"""
struct ExtinctionItpNLTE{
    T,
    ITP_2D <: Interpolations.AbstractInterpolation{T, 2},
    ITP_1D <: Interpolations.AbstractInterpolation{T, 1}
}  <: AbstractBroadening{T}
    σ_atoms::ITP_2D
    σ_hminus::ITP_1D
    σ_h2plus::ITP_1D
    σ_h_ff::ITP_1D
    λ::T
end


#=----------------------------------------------------------------------------
            Functions for computing the background extinction.
----------------------------------------------------------------------------=#

"""
    α_cont(
        itp::ExtinctionItpLTE{T},
        temperature::T,
        electron_density::T,
        hydrogen_density::T,
    )

Calculates continuum extinction in m^-1 using interpolant structure.

# Arguments
- `itp`: Interpolant structure for a given wavelength, hydrogen LTE case
- `temperature`: Temperature in K.
- `electron_density`: Electron number density in m^-3.
- `hydrogen_density`: Number density of hydrogen atoms (all stages), unit m^-3.

# Examples
```julia-repl
julia> npts = 100;
julia> log_temp = LinRange(3.3, 5, npts);
julia> log_ne = LinRange(15, 23, npts);
julia> H_3 = read_atom("H_3.yaml");
julia> itp_lte = create_σ_itp_LTE(500., log_temp, log_ne, H_3,
                                  background_atoms, atom_interpolants);
julia> α_cont(itp_lte, 6000., 1e20, 1e20)
3.59972810095902e-8
```
"""
function α_cont(
    itp::ExtinctionItpLTE{T},
    temperature::T,
    electron_density::T,
    hydrogen_density::T,
)::T where T <: AbstractFloat
    log_temp = log10(temperature)
    log_ne = log10(electron_density)
    α = itp.σ_H(log_temp, log_ne) * hydrogen_density
    α += itp.σ_H2(log_temp, log_ne) * hydrogen_density ^ 2
    α += σ_THOMSON * electron_density
    return α
end


"""
    α_cont(
        itp::ExtinctionItpNLTE{T},
        temperature::T,
        electron_density::T,
        h_neutral_density::T,
        proton_density::T,
    )

Calculates continuum extinction in m^-1 using interpolant structure.

# Arguments
- `itp`: Interpolant structure for a given wavelength, hydrogen NLTE case
- `temperature`: Temperature in K.
- `electron_density`: Electron number density in m^-3.
- `h_neutral_density`: Number density of neutral hydrogen in m^-3.
- `proton_density`: Proton number density in m^-3.

# Examples
```julia-repl
julia> npts = 100;
julia> log_temp = LinRange(3.3, 5, npts);
julia> log_ne = LinRange(15, 23, npts);
julia> itp_nlte = create_σ_itp_LTE(500., log_temp, log_ne,
                                  background_atoms, atom_interpolants);
julia> α_cont(itp_nlte, 6000., 1e20, 1e20, 4.2462e15)
3.5998540603635895e-8
```
"""
function α_cont(
    itp::ExtinctionItpNLTE{T},
    temperature::T,
    electron_density::T,
    h_neutral_density::T,
    proton_density::T,
)::T where T <: AbstractFloat
    log_temp = log10(temperature)
    log_ne = log10(electron_density)
    hydrogen_density = h_neutral_density + proton_density
    α = itp.σ_atoms(log_temp, log_ne) * hydrogen_density
    α += itp.σ_hminus(log_temp) * electron_density * h_neutral_density
    α += itp.σ_h_ff(log_temp) * electron_density * proton_density
    α += itp.σ_h2plus(log_temp) * proton_density * h_neutral_density
    α += ustrip(σ_rayleigh_h(itp.λ * u"nm")) * h_neutral_density
    α += σ_THOMSON * electron_density
    return α
end



"""
    function α_cont_no_itp(
        λ::T,
        temperature::T,
        electron_density::T,
        h_ground_density::T,
        h_neutral_density::T,
        proton_density::T
    ) where T <: AbstractFloat

Calculates continuum extinction without using an interpolation table. Does not include
bound-free processes from background atoms.

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
function α_cont_no_itp(
    λ::T,
    temperature::T,
    electron_density::T,
    h_neutral_density::T,
    proton_density::T
) where T <: AbstractFloat
    λ *= u"nm"
    temperature *= u"K"
    electron_density *= u"m^-3"
    h_neutral_density *= u"m^-3"
    proton_density *= u"m^-3"
    α = α_hminus_ff(λ, temperature, h_neutral_density,  electron_density)
    α += α_hminus_bf(λ, temperature, h_neutral_density, electron_density)
    α += α_hydrogenic_ff(c_0 / λ, temperature, electron_density, proton_density, 1)
    α += α_h2plus_ff(λ, temperature, h_neutral_density, proton_density)
    α += α_h2plus_bf(λ, temperature, h_neutral_density, proton_density)
    α += α_thomson(electron_density)
    α += α_rayleigh_h(λ, h_neutral_density)
    return ustrip(α |> u"m^-1")
end

#=----------------------------------------------------------------------------
            Helper functions to get σ from background atoms.
----------------------------------------------------------------------------=#

"""
    get_atoms_bf_interpolant(atoms::AbstractVector{AtomicModel})

Returns interpolants for bound-free cross section data multiplied with
abundances for each atom.

# Arguments
- `atoms`: A Vector of AtomicModels, with continua.

# Returns
- `tables::Vector{Vector{Interpolations.FilledExtrapolation}}`:
    Interpolation functions.

# Examples
```julia-repl
julia> ATOM_PATH = "/my/atoms/dir/";
julia> atoms = [
    "Al.yaml",
    "C.yaml",
    "Ca.yaml",
    "Fe.yaml",
    "H_6.yaml",
    "He.yaml",
    "KI_fine.yaml",
    "LiI.yaml",
    "Mg.yaml",
    "N.yaml",
    "Na.yaml",
    "NiI.yaml",
    "O.yaml",
    "S.yaml",
    "Si.yaml",
];
julia> background_atoms = Vector{AtomicModel}(undef, length(atoms))
julia> for (index, atom_file) in enumerate(atoms)
           background_atoms[index] = read_atom(join([ATOM_PATH, atom_file]))
       end
julia> atom_interpolants = get_atoms_bf_interpolant(background_atoms);
```
"""
function get_atoms_bf_interpolant(atoms::AbstractVector{AtomicModel})
    tables = [
        Vector{Interpolations.FilledExtrapolation}(
            undef, length(atom.continua)) for atom in atoms
    ]
    for (i, atom) in enumerate(atoms)
        for j in 1:length(atom.continua)
            tables[i][j] = linear_interpolation(
                atom.continua[j].λ,
                atom.continua[j].σ * ABUNDANCES[atom.element],
                extrapolation_bc=0
            )
        end
    end
    return tables
end


"""
    σH_atoms_bf(
        σ_atom_tables::Vector{Vector{Interpolations.FilledExtrapolation}},
        atoms::AbstractVector{AtomicModel},
        λ::T,
        temperature::T,
        electron_density::T
    ) where T <: AbstractFloat

Compute the bound-free cross-sections per hydrogen atom from bf transitions in model atoms.

# Arguments
- `atom_interpolants`: Interpolation functions from get_atoms_bf_interpolant().
- `atoms`: A vector of AtomicModels with continua.
- `λ`: Wavelength in nm.
- `temperature`: Temperature in K.
- `electron_density`: Number density in m^-3.

# Returns
- `σ_λ`: Total cross section per hydrogen atom. Corrected for stimulated emission.
"""
function σH_atoms_bf(
    atom_interpolants::Vector{Vector{Interpolations.FilledExtrapolation}},
    atoms::AbstractVector{AtomicModel},
    λ::T,
    temperature::T,
    electron_density::T
) where T <: AbstractFloat
    σ_λ = 0
    for (i, atom) in enumerate(atoms)
        populations = saha_boltzmann(atom, temperature, electron_density, one(T))
        for j in 1:length(atom.continua)
            σ_j = atom_interpolants[i][j](λ)
            σ_λ += (σ_j * populations[atom.continua[j].lo])
        end
    end
    return σ_λ * (1 - exp(-HC_K / (λ * temperature)))
end


"""
    σH_continuum(λ::T, temperature::T, electron_density::T, ion_frac::T)

Compute continuum cross sections in m^2 per hydrogen atom.

# Arguments
- `λ`: wavelength in nm.
- `temperature`: temperature in K.
- `electron_density`: number density in m^-3.
- `ion_frac` : ionisation fraction.
"""
function σH_continuum(
    λ::T,
    temperature::T,
    electron_density::T,
    ion_frac::T
)::T where T <: AbstractFloat
    λ_u = λ * u"nm"
    ν_u = c_0 / λ_u
    temp = temperature * u"K"
    ne = electron_density * u"m^-3"
    # these will be multiplied by n_HII
    σ = (σ_hminus_bf(λ_u, temp; recipe="wbr") + σ_hminus_ff(λ_u, temp)) * (1 - ion_frac)
    # these will be multiplied by n_HI
    σ += σ_hydrogenic_ff(ν_u, temp, 1) * ion_frac
    σ *= ne
    σ += σ_rayleigh_h(λ_u) * (1 - ion_frac)
    return ustrip(σ |> u"m^2")
end


#=----------------------------------------------------------------------------
            Helper functions for interpolants
----------------------------------------------------------------------------=#

"""
    create_σ_itp_LTE(
        λ::T,
        log_temp::AbstractVector{T},
        log_ne::AbstractVector{T},
        H_atom::AtomicModel,
        background_atoms::AbstractVector{AtomicModel},
        atom_interpolants::Vector{Vector{Interpolations.FilledExtrapolation}},
    )

Create σ interpolant structure for the case of LTE hydrogen populations.

# Arguments
- `λ`: wavelength in nm.
- `log_temp`: sequence of log10 temperature (K) for the table.
- `log_ne`: sequence of log10 electron density (m^-3) for the table.
- `H_atom`: a hydrogen model atom to compute Saha ionisation fractions
- `background_atoms` : sequence of model atoms
- `atom_interpolants` : sequence of bf interpolants corresponding to each background atom
"""
function create_σ_itp_LTE(
    λ::T,
    log_temp::AbstractRange{T},
    log_ne::AbstractRange{T},
    H_atom::AtomicModel,
    background_atoms::AbstractVector{AtomicModel},
    atom_interpolants::Vector{Vector{Interpolations.FilledExtrapolation}},
) where T <: AbstractFloat
    nT = length(log_temp)
    nNe = length(log_ne)
    table_H = Array{T}(undef, nT, nNe)
    table_H2 = Array{T}(undef, nT, nNe)
    λ_u = λ * u"nm"
    for index in CartesianIndices((nT, nNe))
        iT, ie = Tuple(index)
        temp =  10 ^ log_temp[iT]
        ne = 10 ^ log_ne[ie]
        ion_frac = saha_boltzmann(H_atom, temp, ne, 1.0)[end]
        table_H[iT, ie] = σH_continuum(λ, temp, ne, ion_frac)
        table_H[iT, ie] += σH_atoms_bf(atom_interpolants, background_atoms, λ, temp, ne)
        tmp = (σ_h2plus_ff(λ_u, temp * u"K") + σ_h2plus_bf(λ_u, temp * u"K")) |> u"m^5"
        table_H2[iT, ie] = ustrip(tmp) * (1 - ion_frac) * ion_frac
    end
    H_itp = cubic_spline_interpolation((log_temp, log_ne), table_H, extrapolation_bc=Line())
    H2_itp = cubic_spline_interpolation((log_temp, log_ne), table_H2, extrapolation_bc=Line())
    return ExtinctionItpLTE(H_itp, H2_itp, λ)
end


"""
    create_σ_itp_NLTE(
        λ::T,
        log_temp::AbstractVector{T},
        log_ne::AbstractVector{T},
        background_atoms::AbstractVector{AtomicModel},
        atom_interpolants::Vector{Vector{Interpolations.FilledExtrapolation}},
    )

Create σ interpolant structure for the case of LTE hydrogen populations.

# Arguments
- `λ`: wavelength in nm.
- `log_temp`: sequence of log10 temperature (K) for the table.
- `log_ne`: sequence of log10 electron density (m^-3) for the table.
- `background_atoms` : sequence of model atoms
- `atom_interpolants` : sequence of bf interpolants corresponding to each background atom
"""
function create_σ_itp_NLTE(
    λ::T,
    log_temp::AbstractRange{T},
    log_ne::AbstractRange{T},
    background_atoms::AbstractVector{AtomicModel},
    atom_interpolants::Vector{Vector{Interpolations.FilledExtrapolation}},
) where T <: AbstractFloat
    λ_u = λ * u"nm"
    ν_u = c_0 / λ_u
    nT = length(log_temp)
    nNe = length(log_ne)
    table_atoms = Array{T}(undef, nT, nNe)
    table = Array{T}(undef, nT, 3)
    for index in CartesianIndices((nT, nNe))
        iT, ie = Tuple(index)
        temp =  10 ^ log_temp[iT]
        ne = 10 ^ log_ne[ie]
        table_atoms[iT, ie] = σH_atoms_bf(atom_interpolants, background_atoms, λ, temp, ne)
    end
    atoms_itp = cubic_spline_interpolation((log_temp, log_ne), table_atoms, extrapolation_bc=Line())
    for iT in 1:nT
        temp =  (10 ^ log_temp[iT])u"K"
        σ_hminus = (σ_hminus_bf(λ_u, temp; recipe="wbr") + σ_hminus_ff(λ_u, temp)) |> u"m^5"
        table[iT, 1] = ustrip(σ_hminus)
        table[iT, 2] = ustrip((σ_h2plus_ff(λ_u, temp) + σ_h2plus_bf(λ_u, temp)) |> u"m^5")
        table[iT, 3] = ustrip(σ_hydrogenic_ff(ν_u, temp, 1) |> u"m^5")
    end
    hminus_itp = cubic_spline_interpolation((log_temp,), table[:, 1], extrapolation_bc=Line())
    h2plus_itp = cubic_spline_interpolation((log_temp,), table[:, 2], extrapolation_bc=Line())
    h_ff_itp = cubic_spline_interpolation((log_temp,), table[:, 3], extrapolation_bc=Line())
    return ExtinctionItpNLTE(atoms_itp, hminus_itp, h2plus_itp, h_ff_itp, λ)
end
