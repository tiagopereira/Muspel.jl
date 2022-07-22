"""
Reading functions.
"""


"""
Reads atom in YAML format, returns AtomicModel structure
"""
function read_atom(atom_file; FloatT=Float64, IntT=Int)
    data = YAML.load_file(atom_file)
    element = Symbol(data["element"]["symbol"])
    Z = elements[element].number
    if "atomic_mass" in keys(data["element"])
        mass = _assign_unit(data["element"]["atomic_mass"]) |> u"kg"
    else
        mass = elements[element].atomic_mass |> u"kg"
    end
    levels = data["atomic_levels"]
    nlevels = length(levels)
    level_ids = collect(keys(levels))
    # Load and sort levels
    χ = [_assign_unit(level["energy"]) for (_, level) in levels]
    χ = Transparency.wavenumber_to_energy.(χ) .|> u"J"
    idx = sortperm(χ)  # argsort
    χ = χ[idx]
    g = [level["g"] for (_, level) in levels][idx]
    stage = [level["stage"] for (_, level) in levels][idx]
    label = [level["label"] for (_, level) in levels][idx]
    level_ids = level_ids[idx]
    # Continua
    ncont = length(data["radiative_bound_free"])
    continua = Vector{AtomicContinuum}(undef, ncont)
    for (index, cont) in enumerate(data["radiative_bound_free"])
        continua[index] = read_continuum(cont, χ, stage, level_ids;
                                         FloatT=FloatT, IntT=IntT)
    end
    # Lines
    nlines = length(data["radiative_bound_bound"])
    lines = Vector{AtomicLine}(undef, nlines)
    for (index, line) in enumerate(data["radiative_bound_bound"])
        lines[index] = read_line(line, χ, g, stage, level_ids, label, mass;
                                 FloatT=FloatT, IntT=IntT)
    end
    # ...
    # Collisions
    # ...
    # Package things nicely and convert to types
    return AtomicModel{nlevels, FloatT, IntT}(element, nlevels, nlines, ncont, Z,
                                              ustrip(mass), ustrip.(χ), g, stage,
                                              label, lines, continua)
end


"""
Reads continuum transition data in a Dict read from a YAML-formatted atom file.
Needs level energies χ, ionisation stages, and level ids from atom file.
"""
function read_continuum(cont::Dict, χ, stage, level_ids; FloatT=Float64, IntT=Int)
    λedge, up, lo = _read_transition(cont, χ, level_ids)
    # Explicit and hydrogenic cases
    if "cross_section" in keys(cont)
        tmp = reduce(hcat, cont["cross_section"]["value"])
        λ = (tmp[1, :] .* uparse(cont["cross_section"]["unit"][1])) .|> u"nm"
        σ = (tmp[2, :] .* uparse(cont["cross_section"]["unit"][2])) .|> u"m^2"
        idx = sortperm(λ)
        λ = λ[idx]
        σ = σ[idx]
        nλ = length(λ)
    elseif "cross_section_hydrogenic" in keys(cont)
        tmp = cont["cross_section_hydrogenic"]
        nλ = tmp["nλ"]
        λmin = _assign_unit(tmp["λ_min"]) |> u"nm"
        @assert λmin < λedge "Minimum wavelength not shorter than bf edge"
        σ0 = _assign_unit(tmp["σ_peak"])
        λ = LinRange(λmin, λedge, nλ)
        Z_eff = stage[up] - 1
        n_eff = Transparency.n_eff(χ[up], χ[lo], Z_eff)
        σ = σ_hydrogenic_bf_scaled.(σ0, λ, λedge, Z_eff, n_eff) .|> u"m^2"
    else
        error("Photoionisation cross section data missing")
    end
    return AtomicContinuum{nλ, FloatT, IntT}(
        up, lo, nλ, ustrip(λedge), ustrip.(σ), ustrip.(λ))
end


"""
Reads spectral line data in a Dict read from a YAML-formatted atom file.
Needs level energies χ, ionisation stages, labels, level ids and
atomic mass from atom file.
"""
function read_line(line::Dict, χ, g, stage, level_ids, label, mass;
                   FloatT=Float64, IntT=IntT)
    λ0, up, lo = _read_transition(line, χ, level_ids)
    Aul = _assign_unit(line["γ_rad"]) |> u"s^-1"
    Bul = calc_Bji(λ0, Aul) |> u"m^3 / J"
    Blu = g[up] / g[lo] * Bul
    waves = line["wavelengths"]
    f_value = line["f_value"]
    if "data" in keys(waves)
        λ = _assign_unit(waves["data"]) .|> u"nm"
    elseif "type" in keys(waves)
        if waves["type"] == "RH"
            qcore = waves["qcore"]
            qwing = waves["qwing"]
            nλ = waves["nλ"]
            vξ = _assign_unit(waves["vmicro_char"])
            asymm = waves["asymmetric"]
            λ = calc_λline_RH(λ0, nλ, qcore, qwing, vξ; asymm=asymm)
        elseif waves["type"] == "MULTI"
            q0 = waves["q0"]
            qmax = waves["qmax"]
            nλ = waves["nλ"]
            vξ = _assign_unit(waves["qnorm"])
            λ = calc_λline_MULTI(λ0, nλ, q0, qmax, vξ; asymm=false)
        else
            error("Unrecognised wavelength type")
        end
    end
    prof = lowercase(line["type_profile"])
    prd = false
    if prof in ["voigt", "prd"]
        voigt = true
        if prof == "prd"
            prd = true
        end
    elseif prof in ["gaussian", "gauss", "doppler"]
        voigt = false
    else
        error("Unsupported profile type $prof")
    end
    # Energy of the first ionised stage above upper level
    χ∞ = minimum(χ[stage .== stage[up] + 1])
    (vdW_const, vdW_exp) = _read_vdW(line, mass, χ[up], χ[lo], χ∞, stage[up])
    n_vdW = length(vdW_const)
    quad_stark_const = _read_quadratic_stark(line, mass, χ[up], χ[lo], χ∞, stage[up])

    return AtomicLine{n_vdW, FloatT, IntT}(nλ, ustrip(χ[up]), ustrip(χ[lo]), g[up],
                                           g[lo], ustrip(Aul),ustrip(Blu), ustrip(Bul),
                                           ustrip(λ0), f_value, ustrip.(λ), prd, voigt,
                                           label[up], label[lo], vdW_const, vdW_exp,
                                           quad_stark_const)
end


"""
Calculate line wavelengths using recipe from RH.
"""
function calc_λline_RH(λ0::Unitful.Length{T}, nλ, qcore::T, qwing::T, vξ::Unitful.Velocity{T};
                       asymm=true) where T <: Real
    q_to_λ = convert(typeof(λ0), (λ0 * vξ / c_0))
    nhalf = nλ ÷ 2
    if !asymm
        nhalf *= 2
    end
    if qwing <= 2 * qcore
        β = one(T)
    else
        β = qwing / (2 * qcore)
    end
    y = β + sqrt(β * β + (β - 1) * nhalf + 2 - 3 * β)
    b = 2 * log(y) / (nhalf - 1)
    a = qwing / (nhalf - 2 + y*y)
    Δλ = a * ((1:nhalf) .+ exp.(b * ((1:nhalf) .- 1))) .* q_to_λ
    if asymm
        λ = vcat(reverse(λ0 .- Δλ), λ0 .+ Δλ)
    else
        λ = λ0 .+ Δλ
    end
    return λ
end


"""
Calculate line wavelengths using recipe from MULTI.
"""
function calc_λline_MULTI(λ0::Unitful.Length{T}, nλ, q0::T, qmax::T, vξ::Unitful.Velocity{T};
                          asymm=true) where T <: Real
    ν0 = convert(typeof(one(T) * u"Hz"), (c_0 / λ0))
    q = Vector{T}(undef, nλ)

    ten = 10 * one(T)
    al10 = log(ten)
    half = one(T) / 2
    a = ten ^ (q0 + half)
    xmax = log10(a * max(half, qmax - q0 - half))
    if qmax <= q0
        # Linear spacing
        dq = 2 * qmax / (nλ - 1)
        q[1] = -qmax
        for i in 2:nλ
            q[i] = q[i - 1] + dq
        end
    elseif (qmax >= 0) & (q0 >=0)
        if asymm
            dx = 2 * xmax / (nλ - 1)
            for i in 1:nλ
                x = -xmax + (i - 1) * dx
                x10 = ten ^ x
                q[i] = x + (x10 - 1 / x10) / a
            end
        else
            dx = xmax / (nλ - 1)
            for i in 1:nλ
                x = (i - one(T)) * dx
                x10 = ten ^ x
                # Set negative (contrary to MULTI) to ensure consistency
                # with the RH values, which increase from line centre
                q[i] = -(x + (x10 - one(T)) / a)
            end
        end
    end
    ν = ν0 * (1 .+ q * vξ / c_0)
    λ = sort!(convert.(typeof(λ0), c_0 ./ ν))
    return λ
end

#=----------------------------------------------------------------------------
                            Utility functions
----------------------------------------------------------------------------=#
"""
Return a value with a unit given a dictionary read from a YAML file.
"""
_assign_unit(data::Dict) = data["value"] * uparse(data["unit"])


"""
Helper function to parse transition upper and lower levels and compute
wavelength. To be used with fields from the YAML atom format.
"""
function _read_transition(data::Dict, χ, level_ids)
    nlevels = length(χ)
    i = (1:nlevels)[level_ids .== data["transition"][1]][1]
    j = (1:nlevels)[level_ids .== data["transition"][2]][1]
    λ0 = (h * c_0/ abs(χ[i] - χ[j])) |> u"nm"
    up = max(i, j)
    lo = min(i, j)
    return λ0, up, lo
end


"""
Parse type of van der Waals broadening and return the constant
and temperature exponent.
"""
function _read_vdW_single(data::Dict, mass, χup, χlo, χ∞, Z)
    keys_input = lowercase.(keys(data))
    @assert "type" in keys_input "Missing type of van den Waals broadening"
    type = lowercase(data["type"])
    if type == "unsold"
        if "h_coefficient" in keys_input
            h_scaling = data["h_coefficient"]
        else
            h_scaling = 0.0
        end
        if "he_coefficient" in keys_input
            he_scaling = data["he_coefficient"]
        else
            he_scaling = 0.0
        end
        vdw_const = const_unsold(mass, χup, χlo, χ∞, Z;
                                 H_scaling=h_scaling, He_scaling=he_scaling) * u"m^3 / s"
        vdw_exp = 0.3
    elseif type == "abo"  # Barklem
        # Trusting the units
        α = data["α"]["value"]
        σ = data["σ"]["value"]
        vdw_const = const_barklem(mass, α, σ)
        vdw_exp = (1 - α)/2
    elseif type == "deridder_rensbergen"
        α = data["α"]
        β = data["β"]
        h_mass = elements[:H].atomic_mass |> u"kg"
        # Assuming only perturbation by hydrogen
        vdw_const = const_deridder_rensbergen(mass, h_mass, α, β)
        vdw_exp = β
    else
        error("Unsupported van der Waals broadening type: $type")
    end
    return (ustrip(vdw_const |> u"m^3 / s"), vdw_exp)
end


function _read_vdW(line, mass, χup, χlo, χ∞, Z)
    if "broadening_vanderwaals" in keys(line)
        data = line["broadening_vanderwaals"]
        if typeof(data) <: Dict
            data = [data]
        end
        nprocs = length(data)
        arr_const = []
        arr_exp = []
        for process in data
            tmp_const, tmp_exp = _read_vdW_single(process, mass, χup, χlo, χ∞, Z)
            append!(arr_const, tmp_const)
            append!(arr_exp, tmp_exp)
        end
        return (SVector{nprocs, Float64}(arr_const), SVector{nprocs, Float64}(arr_exp))
    else
        return (SVector{0, Float64}(), SVector{0, Float64}())
    end
end


"""
Parse quadratic Stark broadening and return the multiplicative constant.
"""
function _read_quadratic_stark(data::Dict, mass, χup, χlo, χ∞, Z)
    if "broadening_stark" in keys(data)
        coefficient = data["broadening_stark"]["coefficient"]
        if "C_4" in keys(data)  # C_4 provided explicitly
            C_4 = _assign_unit(data["C4"])
        else                    # Use C_4 recipe from Traving 1960
            C_4 = const_quadratic_stark(mass, χup, χlo, χ∞, Z)
        end
        return ustrip((coefficient * C_4) |> u"m^3 / s")
    else
        return 0.0
    end
end
