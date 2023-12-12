var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = Muspel","category":"page"},{"location":"#Muspel","page":"Home","title":"Muspel","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for Muspel.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [Muspel]","category":"page"},{"location":"#Muspel.Atmosphere1D","page":"Home","title":"Muspel.Atmosphere1D","text":"Type for 1D atmospheres. Can contain both 1D only, or 1.5D atmospheres (multiple columns of 1D atmospheres).\n\n\n\n\n\n","category":"type"},{"location":"#Muspel.ExtinctionItpLTE","page":"Home","title":"Muspel.ExtinctionItpLTE","text":"Interpolant structure for continuum extinction, for use when hydrogen populations are computed in LTE (simplifies calculations).\n\n\n\n\n\n","category":"type"},{"location":"#Muspel.ExtinctionItpNLTE","page":"Home","title":"Muspel.ExtinctionItpNLTE","text":"Interpolant structure for continuum extinction, for use when hydrogen populations are given explicitly (e.g. non-equilibrium ionisation or NLTE).\n\n\n\n\n\n","category":"type"},{"location":"#Muspel._assign_unit-Tuple{Dict}","page":"Home","title":"Muspel._assign_unit","text":"Return a value with a unit given a dictionary read from a YAML file.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel._read_broadening_single-Tuple{Dict, Vararg{Any, 5}}","page":"Home","title":"Muspel._read_broadening_single","text":"Reads individual broadening mechanisms and converts each to a multiplier constant and temperature exponent based on different recipes.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel._read_transition-Tuple{Dict, Any, Any}","page":"Home","title":"Muspel._read_transition","text":"Helper function to parse transition upper and lower levels and compute wavelength. To be used with fields from the YAML atom format.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel._spline_coeffs-Tuple{Any, Any}","page":"Home","title":"Muspel._spline_coeffs","text":"Computes coefficients for cubic spline interpolation.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel._spline_stencil-Tuple{Any, Any, Any}","page":"Home","title":"Muspel._spline_stencil","text":"Gets stencil coordinates for cubic spline.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel._w2-Tuple{T} where T<:AbstractFloat","page":"Home","title":"Muspel._w2","text":"Computes weights for linear integration of source function, approximating exp(-Δτ) for very small and very large values of Δτ.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.blackbody_λ-Union{Tuple{T2}, Tuple{T1}, Tuple{T1, T2}} where {T1<:Real, T2<:AbstractFloat}","page":"Home","title":"Muspel.blackbody_λ","text":"Calculates the Blackbody (Planck) function per wavelength in nm and temperature in K. Outputs in kW m^-2 nm^-1.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.calc_line_1D!-Union{Tuple{T}, Tuple{AtomicLine, RTBuffer{T}, Atmosphere1D{1, T, A, V} where {A<:AbstractVector{T}, V<:AbstractVector{T}}, AbstractVector{T}, AbstractVector{T}, ExtinctionItpNLTE{var\"#s37\", T2, ITP_2D, ITP_1D} where {var\"#s37\"<:Real, T2, ITP_2D<:(Interpolations.AbstractInterpolation{var\"#s37\", 2}), ITP_1D<:(Interpolations.AbstractInterpolation{var\"#s37\", 1})}, Interpolations.AbstractInterpolation{<:Real, 2}}} where T<:AbstractFloat","page":"Home","title":"Muspel.calc_line_1D!","text":"function calc_line_1D!(\n    line::AtomicLine,\n    buf::RTBuffer{T},\n    atm::Atmosphere1D{T},\n    n_up::AbstractVector{T},\n    n_lo::AbstractVector{T},\n    σ_itp::ExtinctionItpNLTE{<:Real},\n    voigt_itp::Interpolations.AbstractInterpolation{<:Real, 2},\n)\n\nCalculate emerging disk-centre intensity for a given line in a 1D atmosphere.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.calc_λline_MULTI-Union{Tuple{T}, Tuple{Union{Unitful.Quantity{T, 𝐋, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋, U}} where {L, S}} where U, Any, T, T, Union{Unitful.Quantity{T, 𝐋 𝐓^-1, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋 𝐓^-1, U}} where {L, S}} where U}} where T<:Real","page":"Home","title":"Muspel.calc_λline_MULTI","text":"Calculate line wavelengths using recipe from MULTI.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.calc_λline_RH-Union{Tuple{T}, Tuple{Union{Unitful.Quantity{T, 𝐋, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋, U}} where {L, S}} where U, Any, T, T, Union{Unitful.Quantity{T, 𝐋 𝐓^-1, U}, Unitful.Level{L, S, Unitful.Quantity{T, 𝐋 𝐓^-1, U}} where {L, S}} where U}} where T<:Real","page":"Home","title":"Muspel.calc_λline_RH","text":"Calculate line wavelengths using recipe from RH.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.calc_τ_cont!-Union{Tuple{T}, Tuple{Atmosphere1D{1, T, A, V} where {A<:AbstractVector{T}, V<:AbstractVector{T}}, AbstractVector{<:Real}, ExtinctionItpNLTE{var\"#s39\", T2, ITP_2D, ITP_1D} where {var\"#s39\"<:Real, T2, ITP_2D<:(Interpolations.AbstractInterpolation{var\"#s39\", 2}), ITP_1D<:(Interpolations.AbstractInterpolation{var\"#s39\", 1})}}} where T<:AbstractFloat","page":"Home","title":"Muspel.calc_τ_cont!","text":"Calculate continuum optical depth in the vertical direction, from the observer to the stellar interior. The wavelength is defined by σ_itp.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.create_σ_itp_LTE-Union{Tuple{T}, Tuple{Real, AbstractRange{T}, AbstractRange{T}, AtomicModel, AbstractVector{AtomicModel}, Vector{Vector{Interpolations.FilledExtrapolation}}}} where T<:Real","page":"Home","title":"Muspel.create_σ_itp_LTE","text":"create_σ_itp_LTE(\n    λ::T,\n    log_temp::AbstractVector{T},\n    log_ne::AbstractVector{T},\n    H_atom::AtomicModel,\n    background_atoms::AbstractVector{AtomicModel},\n    atom_interpolants::Vector{Vector{Interpolations.FilledExtrapolation}},\n)\n\nCreate interpolant structure continuum cross section, for the case when hydrogen populations are not given explicitly (and are calculated using Saha). Includes cross sections from bound-free transition present in background atoms, plus the following sources of extinction:\n\n* sources from σH_continuum\n* H ff\n\nArguments\n\nλ: wavelength in nm.\nlog_temp: sequence of log10 temperature (K) for the table.\nlog_ne: sequence of log10 electron density (m^-3) for the table.\nH_atom: a hydrogen model atom to compute Saha ionisation fractions\nbackground_atoms : sequence of model atoms\natom_interpolants : sequence of bf interpolants corresponding to each background atom\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.create_σ_itp_NLTE-Union{Tuple{T}, Tuple{Real, AbstractRange{T}, AbstractRange{T}, AbstractVector{AtomicModel}, Vector{Vector{Interpolations.FilledExtrapolation}}}} where T<:Real","page":"Home","title":"Muspel.create_σ_itp_NLTE","text":"create_σ_itp_NLTE(\n    λ::Real,\n    log_temp::AbstractVector{T},\n    log_ne::AbstractVector{T},\n    background_atoms::AbstractVector{AtomicModel},\n    atom_interpolants::Vector{Vector{Interpolations.FilledExtrapolation}},\n)\n\nCreate interpolant structure continuum cross section, for the case of explicit hydrogen populations. Includes cross sections from bound-free transition present in background atoms, plus the following sources of extinction:\n\nHminus bf and ff\nH2+ molecule bf and ff\nH ff\n\nArguments\n\nλ: wavelength in nm.\nlog_temp: sequence of log10 temperature (K) for the table.\nlog_ne: sequence of log10 electron density (m^-3) for the table.\nbackground_atoms : sequence of model atoms\natom_interpolants : sequence of bf interpolants corresponding to each background atom\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.damping-Union{Tuple{T}, Tuple{T, Real, Real}} where T<:Real","page":"Home","title":"Muspel.damping","text":"Damping constant for γ in rad / s, λ and ΔλD in nm.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.doppler_width-Union{Tuple{T}, Tuple{Any, Any, T}} where T<:AbstractFloat","page":"Home","title":"Muspel.doppler_width","text":"Doppler width for mass in kg, temperature in K\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.feautrier-Union{Tuple{T}, Tuple{AbstractVector{T}, AbstractVector{T}, AbstractVector{T}}} where T<:AbstractFloat","page":"Home","title":"Muspel.feautrier","text":"function feautrier(\n    z::Array{<:Unitful.Length{T}, 1},\n    α::Array{<:PerLength{T}, 1},\n    source_function::Array{<:Unitful.Quantity, 1}\n) where T <: AbstractFloat\n\nCalculate solution to radiative transfer equation using the Feautrier method. Uses algorithm from Rybicki & Hummer, 1991, A&A 245. Returns the height-dependent Feautrier variable P:\n\nP equiv 12 (I^+ + I^-)\n\nCurrently operates under the following assumptions:\n\nThe first index of the variables is the top of the atmosphere\nThe boundary conditions are zero radiation at the top and source function at the bottom\n\nTherefore, the emergent intensity is 2 * P[1], since I^-1=0.\n\nNot properly tested, use with care!\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.get_atoms_bf_interpolant-Tuple{AbstractVector{AtomicModel}}","page":"Home","title":"Muspel.get_atoms_bf_interpolant","text":"get_atoms_bf_interpolant(atoms::AbstractVector{AtomicModel})\n\nReturns interpolants for bound-free cross section data multiplied with abundances for each atom.\n\nArguments\n\natoms: A Vector of AtomicModels, with continua.\n\nReturns\n\ntables::Vector{Vector{Interpolations.FilledExtrapolation}}:   Interpolation functions.\n\nExamples\n\njulia> ATOM_PATH = \"/my/atoms/dir/\";\njulia> atoms = [\n    \"Al.yaml\",\n    \"C.yaml\",\n    \"Ca.yaml\",\n    \"Fe.yaml\",\n    \"H_6.yaml\",\n    \"He.yaml\",\n    \"KI.yaml\",\n    \"Mg.yaml\",\n    \"N.yaml\",\n    \"Na.yaml\",\n    \"NiI.yaml\",\n    \"O.yaml\",\n    \"S.yaml\",\n    \"Si.yaml\",\n];\njulia> background_atoms = Vector{AtomicModel}(undef, length(atoms))\njulia> for (index, atom_file) in enumerate(atoms)\n           background_atoms[index] = read_atom(join([ATOM_PATH, atom_file]))\n       end\njulia> atom_interpolants = get_atoms_bf_interpolant(background_atoms);\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.get_σ_itp-Tuple{Muspel.AbstractAtmos, Real, Vector{String}}","page":"Home","title":"Muspel.get_σ_itp","text":"get_σ_itp(atmos::Atmosphere, λ::Real, atom_files::Vector{String}; npts=100)\n\nConstruct monochromatic continuum cross section interpolant for a given atmosphere, and wavelength λ in nm. Includes all the processes included in create_σ_itp_NLTE, plus all the bound-free cross sections present in all model atom files in the list atom_files.\n\nHere atmos is used only to get the minimum and maximum values of temperature and electron density, to build the interpolation table. The number of points in the table (both for log(temperature) and log(electron density)) is given by npts (default 100).\n\nReturns\n\nExtinctionItpNLTE: continuum cross section interpolant for explicit hydrogen populations. To be used in function α_cont.\n\nExamples\n\njulia> ATOM_PATH = AtomicData.get_atom_dir();\njulia> bckgr_atoms = [\n    \"Al.yaml\",\n    \"C.yaml\",\n    \"Ca.yaml\",\n    \"Fe.yaml\",\n    \"H_6.yaml\",\n    \"He.yaml\",\n    \"KI.yaml\",\n    \"Mg.yaml\",\n    \"N.yaml\",\n    \"Na.yaml\",\n    \"NiI.yaml\",\n    \"O.yaml\",\n    \"S.yaml\",\n    \"Si.yaml\",\n];\njulia> atom_files = [joinpath(ATOM_PATH, a) for a in bckgr_atoms];\njulia> atmos = atmos = read_atmos_rh(MY_ATMOS);\njulia> itp = get_σ_itp(atmos, 500.0, atom_files)\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.h_ionfrac_saha-Union{Tuple{T}, Tuple{T, Real}} where T<:Real","page":"Home","title":"Muspel.h_ionfrac_saha","text":"h_ionfrac_saha(temp::T, electron_density::T)::T where {T <: Real}\n\nCalculate ionisation fraction of hydrogen using Saha.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.incline_atmos-Tuple{Muspel.AbstractAtmos3D, Real, Real}","page":"Home","title":"Muspel.incline_atmos","text":"incline_atmos(atmos_in::AbstractAtmos3D, μ::Real, φ::Real)\n\nTransforms a 3D atmosphere into an inclined coordinate system, given by a rotation by a polar angle θ, given by μ = cos(θ) and an azimuthal angle φ. The output is an atmosphere of the same type, where the quantities have been interpolated to the same number of depth points and the vector quantities projected onto the new axes.\n\nAssumes the following:\n\nOrder of the axes in 3D arrays is (z, y, x)\nHorizontally periodic boundary conditions\nThe first index in the height direction is fixed in the polar rotation\nA right handed system, so that:\n\n      z ↑.            - θ is the clockwise polar rotation from the z axis\n        |  .          - φ is the clockwise azimuthal rotation from the x axis\n        |    .\n        |      .\n        |        *\n        |      . .\n        | θ  .   .\n        |  .     .\n        |.       .\n        +--------.-----→ y\n       /  .      .   /\n      /  φ  .    .\n     /         . . /\n    /- - - - - - .\n  x\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.incline_data!-Tuple{AbstractArray{<:Real, 3}, AbstractVector, Vararg{Real, 4}}","page":"Home","title":"Muspel.incline_data!","text":"incline_data!(\n        data::AbstractArray{<: Real, 3},\n        z::AbstractVector,\n        dx::Real,\n        dy::Real,\n        μ::Real,\n        ϕ::Real\n)\n\nTransforms a 3D array into an inclined coordinate system, according to a polar angle given by μ = cos(θ), and an azimuthal angle ϕ. Uses cubic spline interpolation, and updates data in-place.\n\nThis version only works if height is the first dimension in data.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.incline_data_inv!-Tuple{AbstractArray{<:Real, 3}, AbstractVector, Vararg{Real, 4}}","page":"Home","title":"Muspel.incline_data_inv!","text":"incline_data_inv!(\n        data::AbstractArray{<: Real, 3},\n        z::AbstractVector,\n        dx::Real,\n        dy::Real,\n        μ::Real,\n        ϕ::Real\n)\n\nTransforms a 3D array into an inclined coordinate system, according to a polar angle given by μ = cos(θ), and an azimuthal angle ϕ. Uses cubic spline interpolation, and updates data in-place. Based on trnslt.f90 from Åke Nordlund.\n\nThis version only works if height is the last dimension in data. This is 3-5x faster than incline_data! because the loop order is optimised for this memory configuration.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.piecewise_1D_linear!-Union{Tuple{T}, NTuple{4, AbstractVector{T}}} where T<:AbstractFloat","page":"Home","title":"Muspel.piecewise_1D_linear!","text":"piecewise_1D_linear!(\n    z::AbstractVector{T},\n    α::AbstractVector{T},\n    source_function::AbstractVector{T},\n    intensity::AbstractVector{T};\n    to_end::Bool=false,\n    initial_condition=:source\n) where T <: AbstractFloat\n\nCompute piecewise integration of the radiative transfer equation, assuming linear integration of the source function, for a given height z, extinction α and source_function and intensity (existing array where the output will be saved into). The optional keyword argument to_end defines the direction of the integration: if false (default) will start to integrate intensity from the last element to the first, and if true will integrate from the first element to the last. initial_condition can take two values: :zero for no radiation, or :source (default) to take the source function at the starting point.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.piecewise_1D_linear-Union{Tuple{T}, Tuple{AbstractVector{T}, AbstractVector{T}, AbstractVector{T}}} where T<:AbstractFloat","page":"Home","title":"Muspel.piecewise_1D_linear","text":"piecewise_1D_linear(\n    z::AbstractVector{T},\n    α::AbstractVector{T},\n    source_function::AbstractVector{T};\n    to_end::Bool=false,\n    initial_condition=:source\n) where T <: AbstractFloat\n\nDEPRECATED, use piecewise1Dlinear! Compute piecewise integration of the radiative transfer equation, assuming linear integration of the source function, for a given height z, extinction α and source_function. The optional keyword argument to_end defines the direction of the integration: if false (default) will start to integrate intensity from the last element to the first, and if true will integrate from the first element to the last. initial_condition can take two values: :zero for no radiation, or :source (default) to take the source function at the starting point.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.piecewise_1D_nn-Union{Tuple{T}, Tuple{AbstractVector{T}, AbstractVector{T}, AbstractVector{T}}} where T<:AbstractFloat","page":"Home","title":"Muspel.piecewise_1D_nn","text":"piecewise_1D_nn(\n    z::AbstractVector{T},\n    α::AbstractVector{T},\n    source_function::AbstractVector{T};\n    to_end::Bool=false,\n    initial_condition=:source\n) where T <: AbstractFloat\n\nCompute piecewise integration of the radiative transfer equation, assuming nearest-neighbour integration of the source function, for a given height z, extinction α and source_function. The optional keyword argument to_end defines the direction of the integration: if false (default) will start to integrate intensity from the last element to the first, and if true will integrate from the first element to the last. initial_condition can take two values: :zero for no radiation, or :source (default) to take the source function at the starting point.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.project_vector!-Union{Tuple{A}, Tuple{A, A, A, Real, Real}} where A<:(AbstractArray{<:Real})","page":"Home","title":"Muspel.project_vector!","text":"project_vector!(vx::A, vy::A, vz::A, μ::Real, φ::Real)\n\nProjects a vector in 3D space according to the rotation by an polar angle θ, given by μ = cos(θ) and an azimuthal angle φ. The inputs are the vector components in the x, y, and z axes, for an array of a given size (typically 3D).\n\nAssumes the following:\n\nθ is the clockwise polar rotation from the z axis\nφ is the clockwise azimuthal rotation from the x axis\nThe rotation of the vector is the same as the rotation from the axes (x, y, z) into a new system (x', y', z') given by the rotation matrix:\n\n               Polar          Azimuthal\n   ⎡x'⎤ = ⎡cosθ  0  -sinθ⎤⎡ cosφ  sinφ  0⎤⎡x⎤\n   ⎜y'⎥   ⎜  0   1    0  ⎥⎜-sinφ  cosφ  0⎥⎜y⎥\n   ⎣z'⎦   ⎣sinθ  0   cosθ⎦⎣   0     0   1⎦⎣z⎦\n\n   ⌈x'⎤ = ⎡cosθcosφ  cosθsinφ  -sinφ⎤⎡x⎤\n   |y'⎥   ⎜ -sinφ      cosφ       0 ⎥⎜y⎥\n   ⌊z'⎦   ⎣sinθcosφ  sinθsinφ   cosθ⎦⎣z⎦\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.read_atmos_hpops_multi3d-Tuple{Any, Any, Any}","page":"Home","title":"Muspel.read_atmos_hpops_multi3d","text":"Reads atmosphere in the input format of MULTI3D, at the same time as the hydrogen populations. Only works for a H NLTE run. Only Float32 files are supported at the moment.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.read_atmos_hpops_rh-Tuple{Any, Any}","page":"Home","title":"Muspel.read_atmos_hpops_rh","text":"Reads RH atmosphere. Returns always in single precision.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.read_atmos_multi3d-Tuple{Any, Any}","page":"Home","title":"Muspel.read_atmos_multi3d","text":"Reads atmosphere in the input format of MULTI3D. Only Float32 atmospheres are supported at the moment.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.read_atmos_rh-Tuple{Any}","page":"Home","title":"Muspel.read_atmos_rh","text":"Reads RH atmosphere. Returns always in single precision.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.read_atom-Tuple{Any}","page":"Home","title":"Muspel.read_atom","text":"Reads atom in YAML format, returns AtomicModel structure\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.read_continuum-Tuple{Dict, Any, Any, Any}","page":"Home","title":"Muspel.read_continuum","text":"Reads continuum transition data in a Dict read from a YAML-formatted atom file. Needs level energies χ, ionisation stages, and level ids from atom file.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.read_line-Tuple{Dict, Vararg{Any, 6}}","page":"Home","title":"Muspel.read_line","text":"Reads spectral line data in a Dict read from a YAML-formatted atom file. Needs level energies χ, ionisation stages, labels, level ids and atomic mass from atom file.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.read_mesh-Tuple{Any}","page":"Home","title":"Muspel.read_mesh","text":"Reads mesh file from Bifrost or MULTI3D.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.read_pops_multi3d-NTuple{5, Any}","page":"Home","title":"Muspel.read_pops_multi3d","text":"Reads NLTE populations from MULTI3D output. Does NOT permute dims. Only Float32 files are supported at the moment.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.read_pops_rh-Tuple{Any, Any}","page":"Home","title":"Muspel.read_pops_rh","text":"Reads array with populations for a given species.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.saha_boltzmann!-Union{Tuple{T}, Tuple{N}, Tuple{StaticArraysCore.SVector{N, <:Real}, StaticArraysCore.SVector{N, <:Real}, StaticArraysCore.SVector{N, <:Real}, Real, Real, Real, AbstractVector{T}}} where {N, T<:AbstractFloat}","page":"Home","title":"Muspel.saha_boltzmann!","text":"function saha_boltzmann!(\n    χ::SVector,\n    g::SVector,\n    stage::SVector,\n    temperature::T,\n    electron_density::T,\n    atom_density::T,\n    populations::AbstractArray{T, 1},\n) where T <: AbstractFloat\n\nfunction saha_boltzmann!(\n    atom::AtomicModel,\n    temperature::T,\n    electron_density::T,\n    atom_density::T,\n    populations::AbstractArray{T, 1},\n) where T <: AbstractFloat\n\nInplace version of saha_boltzmann. Calculates atomic level populations according to the Saha-Boltzmann distribution, placing them in an existingpopulations` array.\n\nArguments\n\nχ: level energies in J\ng: statistical weights of levels\nstage: ionisation stage of each level (starting at 1 for neutral)\nOr, instead of the three above, an instance of AtomicModel\ntemperature: temperature in Kelvin\nelectron_density: electron density in m^-3\natom_density: total number density (in all levels) of target species\npopulations: 1D array for output, must be same length as number of levels\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.saha_boltzmann-Union{Tuple{T}, Tuple{N}, Tuple{StaticArraysCore.SVector{N, <:Real}, StaticArraysCore.SVector{N, <:Real}, StaticArraysCore.SVector{N, <:Real}, T, T, Real}} where {N, T<:AbstractFloat}","page":"Home","title":"Muspel.saha_boltzmann","text":"function saha_boltzmann(\n    χ::SVector,\n    g::SVector,\n    stage::SVector,\n    temperature::T,\n    electron_density::T,\n    atom_density::T\n) where T <: AbstractFloat\n\nfunction saha_boltzmann(\n    atom::AtomicModel,\n    temperature::T,\n    electron_density::T,\n    atom_density::T,\n) where T <: AbstractFloat\n\nCalculates atomic level populations according to the Saha-Boltzmann distribution.\n\nArguments\n\nχ: level energies in J\ng: statistical weights of levels\nstage: ionisation stage of each level (starting at 1 for neutral)\nOr, instead of the three above, an instance of AtomicModel\ntemperature: temperature in Kelvin\nelectron_density: electron density in m^-3\natom_density: total number density (in all levels) of target species\n\nReturns\n\npopulations: MVector{nlevels} with relative level populations in m^-3\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.α_cont-Union{Tuple{T}, Tuple{ExtinctionItpLTE{var\"#s42\", T2, ITP_2D} where {var\"#s42\"<:Real, T2, ITP_2D<:(Interpolations.AbstractInterpolation{var\"#s42\", 2})}, T, T, T}} where T<:AbstractFloat","page":"Home","title":"Muspel.α_cont","text":"α_cont(\n    itp::ExtinctionItpLTE{T},\n    temperature::T,\n    electron_density::T,\n    hydrogen_density::T,\n)\n\nCalculates continuum extinction in m^-1 using interpolant structure.\n\nArguments\n\nitp: Interpolant structure for a given wavelength, hydrogen LTE case\ntemperature: Temperature in K.\nelectron_density: Electron number density in m^-3.\nhydrogen_density: Number density of hydrogen atoms (all stages), unit m^-3.\n\nExamples\n\njulia> npts = 100;\njulia> log_temp = LinRange(3.3, 5, npts);\njulia> log_ne = LinRange(15, 23, npts);\njulia> H_3 = read_atom(\"H_3.yaml\");\njulia> itp_lte = create_σ_itp_LTE(500., log_temp, log_ne, H_3,\n                                  background_atoms, atom_interpolants);\njulia> α_cont(itp_lte, 6000., 1e20, 1e20)\n3.59972810095902e-8\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.α_cont-Union{Tuple{T}, Tuple{ExtinctionItpNLTE{var\"#s42\", T2, ITP_2D, ITP_1D} where {var\"#s42\"<:Real, T2, ITP_2D<:(Interpolations.AbstractInterpolation{var\"#s42\", 2}), ITP_1D<:(Interpolations.AbstractInterpolation{var\"#s42\", 1})}, Vararg{T, 4}}} where T<:AbstractFloat","page":"Home","title":"Muspel.α_cont","text":"α_cont(\n    itp::ExtinctionItpNLTE{T},\n    temperature::T,\n    electron_density::T,\n    h_neutral_density::T,\n    proton_density::T,\n)\n\nCalculates continuum extinction in m^-1 using interpolant structure.\n\nArguments\n\nitp: Interpolant structure for a given wavelength, hydrogen NLTE case\ntemperature: Temperature in K.\nelectron_density: Electron number density in m^-3.\nh_neutral_density: Number density of neutral hydrogen in m^-3.\nproton_density: Proton number density in m^-3.\n\nExamples\n\njulia> npts = 100;\njulia> log_temp = LinRange(3.3, 5, npts);\njulia> log_ne = LinRange(15, 23, npts);\njulia> itp_nlte = create_σ_itp_NLTE(500., log_temp, log_ne,\n                                  background_atoms, atom_interpolants);\njulia> α_cont(itp_nlte, 6000., 1e20, 1e20, 4.2462e15)\n3.5998540603635895e-8\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.α_cont_no_itp-Union{Tuple{T}, NTuple{5, T}} where T<:AbstractFloat","page":"Home","title":"Muspel.α_cont_no_itp","text":"function α_cont_no_itp(\n    λ::T,\n    temperature::T,\n    electron_density::T,\n    h_ground_density::T,\n    h_neutral_density::T,\n    proton_density::T\n) where T <: AbstractFloat\n\nCalculates continuum extinction without using an interpolation table. Does not include bound-free processes from background atoms.\n\nArguments\n\nλ: Wavelength in nm.\ntemperature: Temperature in K.\nelectron_density: Electron number density in m^-3.\nh_ground_density: Number density of hydrogen in the ground state, unit m^-3.\nh_neutral_density: Number density of neutral hydrogen, unit m^-3.\nproton_density: Proton number density in m^-3.\n\nReturns\n\nα: Continuous extinction (Float) in m^-1.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.σH_atoms_bf-Union{Tuple{T}, Tuple{Vector{Vector{Interpolations.FilledExtrapolation}}, AbstractVector{AtomicModel}, Real, T, T}} where T<:AbstractFloat","page":"Home","title":"Muspel.σH_atoms_bf","text":"σH_atoms_bf(\n    σ_atom_tables::Vector{Vector{Interpolations.FilledExtrapolation}},\n    atoms::AbstractVector{AtomicModel},\n    λ::T,\n    temperature::T,\n    electron_density::T\n) where T <: AbstractFloat\n\nCompute the bound-free cross-sections per hydrogen atom from bf transitions in model atoms.\n\nArguments\n\natom_interpolants: Interpolation functions from getatomsbf_interpolant().\natoms: A vector of AtomicModels with continua.\nλ: Wavelength in nm.\ntemperature: Temperature in K.\nelectron_density: Number density in m^-3.\n\nReturns\n\nσ_λ: Total cross section per hydrogen atom. Corrected for stimulated emission.\n\n\n\n\n\n","category":"method"},{"location":"#Muspel.σH_continuum-Union{Tuple{T}, Tuple{Real, T, T, Real}} where T<:AbstractFloat","page":"Home","title":"Muspel.σH_continuum","text":"σH_continuum(λ::T, temperature::T, electron_density::T, ion_frac::T)\n\nCompute continuum cross sections in m^2 per hydrogen atom.\n\nArguments\n\nλ: wavelength in nm.\ntemperature: temperature in K.\nelectron_density: number density in m^-3.\nion_frac : ionisation fraction.\n\n\n\n\n\n","category":"method"}]
}
