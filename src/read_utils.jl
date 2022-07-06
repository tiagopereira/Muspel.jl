"""
Reading functions.
"""

"""
Return a value with a unit given a dictionary read from a YAML file.
"""
parse_unit(data::Dict) = data["value"] * uparse(data["unit"])


"""
Reads atom in YAML format, returns AtomicModel structure
"""
function read_atom(atom_file)
    data = YAML.load_file(atom_file)
    element = Symbol(data["element"]["symbol"])
    Z = elements[element].number
    mass = ustrip(elements[element].atomic_mass |> u"kg")
    levels = data["atomic_levels"]
    nlevels = length(levels)
    level_ids = collect(keys(levels))
    # Load and sort levels
    χ = [parse_unit(level["energy"]) for (_, level) in levels]
    χ = ustrip(Transparency.wavenumber_to_energy.(χ) .|> u"J")
    idx = sortperm(χ)  # argsort
    FloatT = eltype(χ)
    IntT = typeof(Z)
    χ = SVector{nlevels, FloatT}(χ[idx])
    g = SVector{nlevels, IntT}([level["g"] for (_, level) in levels][idx])
    stage = SVector{nlevels, IntT}([level["stage"] for (_, level) in levels][idx])
    label = [level["label"] for (_, level) in levels][idx]
    level_ids = level_ids[idx]
    # Continua
    ncont = length(data["radiative_bound_free"])
    continua = Vector{AtomicContinuum}(undef, ncont)
    for (index, cont) in enumerate(data["radiative_bound_free"])
        # indices of levels, find upper and lower
        i = (1:nlevels)[level_ids .== cont["transition"][1]][1]
        j = (1:nlevels)[level_ids .== cont["transition"][2]][1]
        λedge = ustrip(h * c_0/ (abs(χ[i] - χ[j]))u"J" |> u"nm")
        up = max(i, j)
        lo = min(i, j)
        # Explicit and hydrogenic cases
        if "cross_section" in keys(cont)
            tmp = reduce(hcat, cont["cross_section"]["value"])
            λ = ustrip.((tmp[1, :] .* uparse(cont["cross_section"]["unit"][1])) .|> u"nm")
            σ = ustrip.((tmp[2, :] .* uparse(cont["cross_section"]["unit"][2])) .|> u"m^2")
            idx = sortperm(λ)
            λ = λ[idx]
            σ = σ[idx]
        elseif "cross_section_hydrogenic" in keys(cont)
            tmp = cont["cross_section_hydrogenic"]
            nλ = tmp["nλ"]
            λmin = ustrip(parse_unit(tmp["λ_min"]) |> u"nm")
            @assert λmin < λedge "Minimum wavelength not shorter than bf edge"
            σ0 = parse_unit(tmp["σ_peak"])
            λ = LinRange(λmin, λedge, nλ)
            Z_eff = stage[up] - 1
            n_eff = Transparency.n_eff(χ[up] * u"J", χ[lo] * u"J", Z_eff)
            σ = hydrogenic_bf_σ_scaled.(σ0, λ * u"nm", λedge * u"nm", Z_eff, n_eff)
            σ = ustrip(σ .|> u"m^2")
        else
            error("Photoionisation cross section data missing")
        end
        continua[index] = AtomicContinuum(λedge, up, lo, σ, λ)
    end
    # Lines
    # ...
    # Collisions
    # ...
    return AtomicModel{nlevels, FloatT, IntT}(element, nlevels, 1, ncont, Z, mass,
                                              χ, g, stage, label, continua)
end
