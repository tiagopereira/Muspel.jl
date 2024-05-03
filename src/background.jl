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
    T1,
    T2,
    ITP_2D <: Interpolations.AbstractInterpolation{T1, 2}
}  <: AbstractExtinctionItp{T1}
    σ_H::ITP_2D
    σ_H2::ITP_2D
    λ::T2
end


"""
Interpolant structure for continuum extinction, for use when hydrogen
populations are given explicitly (e.g. non-equilibrium ionisation or NLTE).
"""
struct ExtinctionItpNLTE{
    T1,
    T2,
    ITP_2D <: Interpolations.AbstractInterpolation{T1, 2},
    ITP_1D <: Interpolations.AbstractInterpolation{T1, 1}
}  <: AbstractExtinctionItp{T1}
    σ_atoms::ITP_2D
    σ_hminus::ITP_1D
    σ_h2plus::ITP_1D
    σ_h_ff::ITP_1D
    λ::T2
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
    itp::ExtinctionItpLTE{<: Real},
    temperature::T,
    electron_density::T,
    hydrogen_density::T,
)::T where T <: AbstractFloat
    log_temp = log10(temperature)
    log_ne = log10(electron_density)
    α = itp.σ_H(log_temp, log_ne) * hydrogen_density
    α += (itp.σ_H2(log_temp, log_ne) * hydrogen_density) * hydrogen_density
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
julia> itp_nlte = create_σ_itp_NLTE(500., log_temp, log_ne,
                                  background_atoms, atom_interpolants);
julia> α_cont(itp_nlte, 6000., 1e20, 1e20, 4.2462e15)
3.5998540603635895e-8
```
"""
function α_cont(
    itp::ExtinctionItpNLTE{<: Real},
    temperature::T,
    electron_density::T,
    h_neutral_density::T,
    proton_density::T,
)::T where T <: AbstractFloat
    log_temp = log10(temperature)
    log_ne = log10(electron_density)
    hydrogen_density = h_neutral_density + proton_density
    α = itp.σ_atoms(log_temp, log_ne) * hydrogen_density
    α += (itp.σ_hminus(log_temp) * electron_density) * h_neutral_density
    α += (itp.σ_h_ff(log_temp) * electron_density) * proton_density
    α += (itp.σ_h2plus(log_temp) * proton_density) * h_neutral_density
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
            Functions to get bf cross sections from background atoms.
----------------------------------------------------------------------------=#

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
    λ::Real,
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
    λ::Real,
    temperature::T,
    electron_density::T,
    ion_frac::Real
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
        Functions to create interpolants of continuum cross section
----------------------------------------------------------------------------=#
"""
    get_σ_itp(atmos::Atmosphere, λ::Real, atom_files::Vector{String}; npts=100)

Construct monochromatic continuum cross section interpolant for a given atmosphere,
and wavelength λ in nm. Includes all the processes included in `create_σ_itp_NLTE`,
plus all the bound-free cross sections present in all model atom files in the
list `atom_files`.

Here `atmos` is used only to get the minimum and maximum values of temperature
and electron density, to build the interpolation table. The number of points in the
table (both for log(temperature) and log(electron density)) is given by `npts`
(default 100).

# Returns
- ExtinctionItpNLTE: continuum cross section interpolant for explicit hydrogen
  populations. To be used in function `α_cont`.

# Examples
```julia-repl
julia> ATOM_PATH = AtomicData.get_atom_dir();
julia> bckgr_atoms = [
    "Al.yaml",
    "C.yaml",
    "Ca.yaml",
    "Fe.yaml",
    "H_6.yaml",
    "He.yaml",
    "KI.yaml",
    "Mg.yaml",
    "N.yaml",
    "Na.yaml",
    "NiI.yaml",
    "O.yaml",
    "S.yaml",
    "Si.yaml",
];
julia> atom_files = [joinpath(ATOM_PATH, a) for a in bckgr_atoms];
julia> atmos = atmos = read_atmos_rh(MY_ATMOS);
julia> itp = get_σ_itp(atmos, 500.0, atom_files)
```
"""
function get_σ_itp(atmos::AbstractAtmos, λ::Real, atom_files::Vector{String}; npts=100)
    bckgr_atoms = Vector{AtomicModel}(undef, length(atom_files))
    for (index, file) in enumerate(atom_files)
        bckgr_atoms[index] = read_atom(file)
    end
    atom_interpolants = get_atoms_bf_interpolant(bckgr_atoms)
    t_range = log10.([Float64(minimum(atmos.temperature)),
                      Float64(maximum(atmos.temperature))])
    ne_range = log10.([Float64(minimum(atmos.electron_density)),
                       Float64(maximum(atmos.electron_density))])
    log_temp = LinRange(minimum(t_range), maximum(t_range), npts)
    log_ne = LinRange(minimum(ne_range), maximum(ne_range), npts)
    return create_σ_itp_NLTE(λ, log_temp, log_ne, bckgr_atoms, atom_interpolants)
end


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
    "KI.yaml",
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
    create_σ_itp_LTE(
        λ::T,
        log_temp::AbstractVector{T},
        log_ne::AbstractVector{T},
        H_atom::AtomicModel,
        background_atoms::AbstractVector{AtomicModel},
        atom_interpolants::Vector{Vector{Interpolations.FilledExtrapolation}},
    )

Create interpolant structure continuum cross section, for the case when hydrogen
populations are not given explicitly (and are calculated using Saha). Includes cross
sections from bound-free transition present in background atoms, plus the
following sources of extinction:

    * sources from σH_continuum
    * H ff

# Arguments
- `λ`: wavelength in nm.
- `log_temp`: sequence of log10 temperature (K) for the table.
- `log_ne`: sequence of log10 electron density (m^-3) for the table.
- `H_atom`: a hydrogen model atom to compute Saha ionisation fractions
- `background_atoms` : sequence of model atoms
- `atom_interpolants` : sequence of bf interpolants corresponding to each background atom
"""
function create_σ_itp_LTE(
    λ::Real,
    log_temp::AbstractRange{T},
    log_ne::AbstractRange{T},
    H_atom::AtomicModel,
    background_atoms::AbstractVector{AtomicModel},
    atom_interpolants::Vector{Vector{Interpolations.FilledExtrapolation}},
) where T <: Real
    nT = length(log_temp)
    nNe = length(log_ne)
    # Tables must be Float64 to avoid loss in accuracy
    table_H = Array{Float64}(undef, nT, nNe)
    table_H2 = Array{Float64}(undef, nT, nNe)
    λ_u = λ * u"nm"
    for index in CartesianIndices((nT, nNe))
        iT, ie = Tuple(index)
        temp = 10 ^ log_temp[iT]
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
        λ::Real,
        log_temp::AbstractVector{T},
        log_ne::AbstractVector{T},
        background_atoms::AbstractVector{AtomicModel},
        atom_interpolants::Vector{Vector{Interpolations.FilledExtrapolation}},
    )

Create interpolant structure continuum cross section, for the case of explicit
hydrogen populations. Includes cross sections from bound-free transition present in
background atoms, plus the following sources of extinction:

* Hminus bf and ff
* H2+ molecule bf and ff
* H ff

# Arguments
- `λ`: wavelength in nm.
- `log_temp`: sequence of log10 temperature (K) for the table.
- `log_ne`: sequence of log10 electron density (m^-3) for the table.
- `background_atoms` : sequence of model atoms
- `atom_interpolants` : sequence of bf interpolants corresponding to each background atom
"""
function create_σ_itp_NLTE(
    λ::Real,
    log_temp::AbstractRange{T},
    log_ne::AbstractRange{T},
    background_atoms::AbstractVector{AtomicModel},
    atom_interpolants::Vector{Vector{Interpolations.FilledExtrapolation}},
) where T <: Real
    λ_u = λ * u"nm"
    ν_u = c_0 / λ_u
    nT = length(log_temp)
    nNe = length(log_ne)
    # Tables must be Float64 to avoid loss in accuracy
    table_atoms = Array{Float64}(undef, nT, nNe)
    table = Array{Float64}(undef, nT, 3)
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
