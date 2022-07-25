var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = Muspel","category":"page"},{"location":"#Muspel","page":"Home","title":"Muspel","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for Muspel.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Muspel]","category":"page"},{"location":"#Muspel.Tables_σ","page":"Home","title":"Muspel.Tables_σ","text":"function Tables_σ(\n    λ::Vector{T},\n    log_temp::StepRangeLen,\n    log_ne::StepRangeLen,\n    atoms::Vector{AtomicModel}\n) where T<:AbstractFloat\n\nStructure for the cross-section interpolation tables used for the background extinction. Each field is an array of CubicSplineInterpolation functions. There is one function per wavelength.\n\nArguments:\n\nλ: Wavelengths for which tables are computed.\nlog_temp: Range of log10 temperatures for the tables.\nlog_ne: Range of log10 electron densities for the tables.\natoms: Background atoms with continua to include.\n\nFields:\n\ntable_nh::Vector{Interpolations.Extrapolation}:           Takes arguments (log10 temperature, log10 electron density). Contains           contribution from atomic continua, H- ff and H- bf. Multiply by hydrogen           density to get extinction per m^-1.\ntable_ne::Vector{Interpolations.Extrapolation}:           Takes argument (log10 temperature). Contains hydrogenic ff contribution from           protons. Multiply by proton and electron density to get extinction per m^-1.\n\n\n\n\n\n","category":"type"},{"location":"#Muspel._assign_unit-Tuple{Dict}","page":"Home","title":"Muspel._assign_unit","text":"Return a value with a unit given a dictionary read from a YAML file.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel._read_quadratic_stark-Tuple{Dict, Any, Any, Any, Any, Any}","page":"Home","title":"Muspel._read_quadratic_stark","text":"Parse quadratic Stark broadening and return the multiplicative constant.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel._read_transition-Tuple{Dict, Any, Any}","page":"Home","title":"Muspel._read_transition","text":"Helper function to parse transition upper and lower levels and compute wavelength. To be used with fields from the YAML atom format.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel._read_vdW_single-Tuple{Dict, Any, Any, Any, Any, Any}","page":"Home","title":"Muspel._read_vdW_single","text":"Parse type of van der Waals broadening and return the constant and temperature exponent.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.calc_λline_MULTI-Union{Tuple{T}, Tuple{Union{Unitful.Quantity{T, 𝐋, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋, U}} where {L, S}} where U, Any, T, T, Union{Unitful.Quantity{T, 𝐋 𝐓^-1, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋 𝐓^-1, U}} where {L, S}} where U}} where T<:Real","page":"Home","title":"Muspel.calc_λline_MULTI","text":"Calculate line wavelengths using recipe from MULTI.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.calc_λline_RH-Union{Tuple{T}, Tuple{Union{Unitful.Quantity{T, 𝐋, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋, U}} where {L, S}} where U, Any, T, T, Union{Unitful.Quantity{T, 𝐋 𝐓^-1, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋 𝐓^-1, U}} where {L, S}} where U}} where T<:Real","page":"Home","title":"Muspel.calc_λline_RH","text":"Calculate line wavelengths using recipe from RH.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.read_atom-Tuple{Any}","page":"Home","title":"Muspel.read_atom","text":"Reads atom in YAML format, returns AtomicModel structure\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.read_continuum-Tuple{Dict, Any, Any, Any}","page":"Home","title":"Muspel.read_continuum","text":"Reads continuum transition data in a Dict read from a YAML-formatted atom file. Needs level energies χ, ionisation stages, and level ids from atom file.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.read_line-Tuple{Dict, Any, Any, Any, Any, Any, Any}","page":"Home","title":"Muspel.read_line","text":"Reads spectral line data in a Dict read from a YAML-formatted atom file. Needs level energies χ, ionisation stages, labels, level ids and atomic mass from atom file.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.saha_boltzmann-Union{Tuple{T}, Tuple{StaticArraysCore.SVector, StaticArraysCore.SVector, StaticArraysCore.SVector, T, T}} where T<:AbstractFloat","page":"Home","title":"Muspel.saha_boltzmann","text":"function saha_boltzmann(χ::SVector,\n                        g::SVector,\n                        stage::SVector,\n                        temperature::T,\n                        electron_density::T) where T <: AbstractFloat\n\nfunction saha_boltzmann(atom::AtomicModel,\n                        temperature::T,\n                        electron_density::T) where T <: AbstractFloat\n\nCalculates atomic level populations according to the Saha-Boltzmann distribution.\n\nArguments\n\nχ: level energies in J\ng: statistical weights of levels\nstage: ionisation stage of each level (starting at 1 for neutral)\nOr, instead of the three above, an instance of AtomicModel\ntemperature: temperature in Kelvin\nelectron_density: electron density in m^-3\n\nReturns\n\npopulations: MVector{nlevels} with relative level populations in m^-3\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.α_atoms_bf-Union{Tuple{T}, Tuple{Vector{Vector{Interpolations.FilledExtrapolation}}, Vector{AtomicModel}, T, T, T, T}} where T<:AbstractFloat","page":"Home","title":"Muspel.α_atoms_bf","text":"function α_atoms_bf(\n    σ_atom_tables::Vector{Vector{Interpolations.FilledExtrapolation}},\n    atoms::Vector{AtomicModel},\n    λ::T,\n    temperature::T,\n    electron_density::T,\n    hydrogen_density::T\n) where T <: AbstractFloat\n\nComputes the extinction per meter from continua in atoms.\n\nArguments\n\nσ_atom_tables: Interpolation tables from σatomsbf_tables().\natoms: Vector of AtomicModels with continua.\nλ : Wavelength in nm.\ntemperature: temperature in K.\nelectron_density: Electron density in m^-3.\nhydrogen_density: Total hydrogen density in m^-3.\n\nReturns\n\nα_λ: Extinction per meter from bound-free transitions (Float).\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.α_cont-Tuple{Union{Unitful.Quantity{T, 𝐋, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋, U}} where {L, S}} where {T, U}, Union{Unitful.Quantity{T, 𝚯, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝚯, U}} where {L, S}} where {T, U}, Union{Unitful.Quantity{T, 𝐋^-3, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋^-3, U}} where {L, S}} where {T, U}, Union{Unitful.Quantity{T, 𝐋^-3, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋^-3, U}} where {L, S}} where {T, U}, Union{Unitful.Quantity{T, 𝐋^-3, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋^-3, U}} where {L, S}} where {T, U}, Union{Unitful.Quantity{T, 𝐋^-3, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋^-3, U}} where {L, S}} where {T, U}}","page":"Home","title":"Muspel.α_cont","text":"function α_cont(\n    λ::Unitful.Length,\n    temperature::Unitful.Temperature,\n    electron_density::NumberDensity,\n    h_ground_density::NumberDensity,\n    h_neutral_density::NumberDensity,\n    proton_density::NumberDensity\n)\n\nfunction α_cont(\n    λ::T,\n    temperature::T,\n    electron_density::T,\n    h_ground_density::T,\n    h_neutral_density::T,\n    proton_density::T\n) where T <: AbstractFloat\n\nfunction α_cont(\n    atoms::Vector{AtomicModel},\n    σ_atom_tables::Vector{Vector{Interpolations.FilledExtrapolation}},\n    λ::T,\n    temperature::T,\n    electron_density::T,\n    h_ground_density::T,\n    h_neutral_density::T,\n    proton_density::T\n) where T <: AbstractFloat\n\nCalculates continuum extinction according to recipe in RH using functions from Transparecy. Optionally includes continua from atom files. Atoms are treated in LTE.\n\nArguments\n\natoms: (optional) Atoms with continua to include.\nσ_atom_tables: (optional) Interpolation functions from σatomsbf_tables().\nλ: Wavelength in nm or with a Unitful.Length unit.\ntemperature: Unit K.\nelectron_density: Number density, unit m^-3.\nh_ground_density: Number density of hydrogen in the ground state.\nh_neutral_density: Number density of neutral hydrogen.\nproton_density: Number density, unit m^-3.\n\nReturns\n\nα: Continuous extinction (Float) in m^-1. Has units if input had units.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.α_cont-Union{Tuple{T}, Tuple{Tables_σ, Integer, T, T, T, T, T, T}} where T<:AbstractFloat","page":"Home","title":"Muspel.α_cont","text":"function α_cont(\n    tables::Tables_σ,\n    iλ::Integer,\n    λ::T,\n    temperature::T,\n    electron_density::T,\n    h_ground_density::T,\n    h_neutral_density::T,\n    proton_density::T\n) where T<: AbstractFloat\n\nfunction α_cont(tables::Tables_σ, iλ::Integer, args...) where args <: Unitful.Quantity\n\nIf given a Tables_σ struct, the function returns the extinction computed from interpolation tables of the cross-sections. This is significally faster than the non-tabulated version when the computation includes many atomic continua.\n\nArguments\n\n- `tables`: Pre-computed interpolation functions for \"cross-sections\".\n- `iλ`: The index in `tables` for the wavelength λ.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.σ_atoms_bf-Union{Tuple{T}, Tuple{Vector{Vector{Interpolations.FilledExtrapolation}}, Vector{AtomicModel}, T, T, T}} where T<:AbstractFloat","page":"Home","title":"Muspel.σ_atoms_bf","text":"function σ_atoms_bf(\n    σ_atom_tables::Vector{Vector{Interpolations.FilledExtrapolation}},\n    atoms::Vector{AtomicModel},\n    λ::T,\n    temperature::T,\n    electron_density::T\n) where T <: AbstractFloat\n\nComputes the bound-free cross-sections from atom files multiplied with the population of the lower level of each bound-free transition relative to the total hydrogen population. To get total bound-free extinction multiply with hydrogen_density.\n\nArguments\n\nσ_atom_tables: Interpolation functions from σatomsbf_tables().\natoms: A vector of AtomicModels with continua.\nλ: Wavelength in nm.\ntemperature: Temperature in K.\nelectron_density: Number density in m^-3.\n\nReturns\n\nσ_λ: Sum of all cross-sections multiplied with relative populations, abundancies and   correction for stimulated emission (Float).\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.σ_atoms_bf_tables-Tuple{Vector{AtomicModel}}","page":"Home","title":"Muspel.σ_atoms_bf_tables","text":"function σ_atoms_bf_tables(atoms::Vector{AtomicModel})\n\nReturns interpolation functions for the bound-free cross section data multiplied with abundances for each atom.\n\nArguments\n\natoms: A Vector of AtomicModels with continua.\n\nReturns\n\nσ_atom_tables::Vector{Vector{Interpolations.FilledExtrapolation}}:   Interpolation functions.\n\n\n\n\n\n","category":"method"}]
}
