"""
Reading functions.
"""


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
    χ = [level["energy"]["value"] for (_, level) in levels]
    idx = sortperm(χ)  # argsort
    unit = uparse.([level["energy"]["unit"] for (_, level) in levels])[idx]
    FloatT = eltype(χ)
    IntT = typeof(Z)
    # Convert to energy, but without units
    χ = SVector{nlevels, FloatT}(
        ustrip(Transparency.wavenumber_to_energy.(χ[idx] .* unit) .|> u"J"))
    g = SVector{nlevels, IntT}([level["g"] for (_, level) in levels][idx])
    stage = SVector{nlevels, IntT}([level["stage"] for (_, level) in levels][idx])
    label = [level["label"] for (_, level) in levels][idx]
    level_ids = level_ids[idx]
    indices = 1:nlevels
    # Continua
    ncont = length(data["radiative_bound_free"])
    continua = Vector{AtomicContinuum}(undef, ncont)
    counter = 1
    for cont in data["radiative_bound_free"]
        tmp = reduce(hcat, cont["cross_section"]["value"])
        λ = ustrip.((tmp[1, :] .* uparse(cont["cross_section"]["unit"][1])) .|> u"nm")
        σ = ustrip.((tmp[2, :] .* uparse(cont["cross_section"]["unit"][2])) .|> u"m^2")
        idx = sortperm(λ)
        λ = λ[idx]
        σ = σ[idx]
        # indices of levels, find upper and lower
        i = (1:nlevels)[level_ids .== cont["transition"][1]][1]
        j = (1:nlevels)[level_ids .== cont["transition"][2]][1]
        up = max(i, j)
        lo = min(i, j)
        continua[counter] = AtomicContinuum(up, lo; kind="explicit", σ=σ, λ=λ)
        counter += 1
    end
    # Lines
    # ...
    # Collisions
    # ...
    AtomicModel{nlevels, FloatT, IntT}(element, nlevels, 1, ncont, Z, mass,
                                       χ, g, stage, label, continua)
end
