"""
Calculates the effective Landé factor.
"""
function g_eff(g1, g2, J1, J2)
    d = J1 * (J1 + 1) - J2 * (J2 + 1)
    return (g1 + g2) * 0.5f0 + (g1 - g2) * d * 0.25f0
end


"""
Calculates Landé factor from LS coupling.
"""
function g_LS(S, L, J)
    if J == 0
        return 1f0
    else
        return 1f0 + 0.5f0 * (J * (J + 1) + S * (S + 1) - L * (L + 1)) / (J * (J + 1))
    end
end


"""
Calculates strengths of Zeeman components. These strengths are normalised.
"""
function zeeman_strength(J_l, J_u, M_l, M_u)
    J, M = J_l, M_l
    if J_l == J_u == 0
        error("Invalid transition, J_u = J_l = 0")
    end
    if M_u - M_l == 1
        if J_u - J_l == 1
            strength = (3 * (J + M + 1) * (J + M + 2)) / (2 * (J + 1) * (2 * J + 1) * (2 * J + 3))
        elseif J_u - J_l == 0
            strength = (3 * (J - M) * (J + M + 1)) / (2 * J * (J + 1) * (2 * J + 1))
        elseif J_u - J_l == -1
            strength = (3 * (J - M) * (J - M - 1)) / (2 * J * (2 * J - 1) * (2 * J + 1))
        else
            error("Invalid transition, J_u - J_l != -1, 0, 1")
        end
    elseif M_u - M_l == 0
        if J_u - J_l == 1
            strength = (3 * (J - M + 1) * (J + M + 1)) / ((J + 1) * (2 * J + 1) * (2 * J + 3))
        elseif J_u - J_l == 0
            strength = 3 * M^2 / (J * (J + 1) * (2 * J + 1))
        elseif J_u - J_l == -1
            strength = (3 * (J - M) * (J + M)) / (J * (2 * J - 1) * (2 * J + 1))
        else
            error("Invalid transition, J_u - J_l != -1, 0, 1")
        end
    elseif M_u - M_l == -1
        if J_u - J_l == 1
            strength = (3 * (J - M + 1) * (J - M + 2)) / (2 * (J + 1) * (2 * J + 1) * (2 * J + 3))
        elseif J_u - J_l == 0
            strength = (3 * (J + M) * (J - M + 1)) / (2 * J * (J + 1) * (2 * J + 1))
        elseif J_u - J_l == -1
            strength = (3 * (J + M) * (J + M - 1)) / (2 * J * (2 * J - 1) * (2 * J + 1))
        else
            error("Invalid transition, J_u - J_l != -1, 0, 1")
        end
    else
        error("Invalid transition, M_u - M_l != -1, 0, 1")
    end
    return strength
end


"""
Parses an atomic level label in RH format, to extract the atomic numbers S and L.
"""
function parse_label_LS(label)
    orbital = Dict(
        'S' => 0, 'P' => 1, 'D' => 2, 'F' => 3, 'G' => 4, 'H' => 5, 'I' => 6,
        'K' => 7, 'L' => 8, 'M' => 9, 'N' => 10, 'O' => 11, 'Q' => 12, 'R' => 13,
        'T' => 14, 'U' => 15,'V' => 11, 'W' => 12, 'X' => 19, 'Z' => 20
    )
    # Remove characters after O and E, handle all strange cases
    if count('E', label) > 0 && count('O', label) > 0
        if findlast("E", label) > findlast("O", label)
            label = split(label, "E")[end-1]
        else
            label = split(label, "O")[end-1]
        end
    elseif count('E', label) > 0
        label = split(label, "E")[end-1]
    elseif count('O', label) > 0
        label = split(label, "O")[end-1]
    end
    tmp = split(label)
    if length(tmp) == 0
        error("Found empty label")
    else
        term = tmp[end]
        try
            multiplicity = parse(Int, term[end-1])
            S = (multiplicity - 1) // 2
            L = orbital[uppercase(term[end])]
            return (S, L)
        catch e
            error("Could not parse term symbol: $term")
        end
    end
end


"""
Computes Zeeman strengths and shifts for the π and σ components.
"""
function get_zeeman_components(Sl, Ll, Jl, Su, Lu, Ju)
    σr_strength = zeros(Float32, 0)
    σr_shift = zeros(Float32, 0)
    π_strength = zeros(Float32, 0)
    π_shift = zeros(Float32, 0)
    σb_strength = zeros(Float32, 0)
    σb_shift = zeros(Float32, 0)

    gu = g_LS(Su, Lu, Ju)
    gl = g_LS(Sl, Ll, Jl)

    # consider only valid cases where ΔJ = +1, 0, -1
    if abs(Jl - Ju) in [-1, 0, 1]
        for Ml in -Jl:Jl, Mu in -Ju:Ju
            ΔM = Mu - Ml
            if ΔM == -1  # sigma^- or r components
                strength = zeeman_strength(Jl, Ju, Ml, Mu)
                if strength > 0
                    push!(σr_strength, strength)
                    push!(σr_shift, gu*Mu - gl*Ml)
                end
            elseif ΔM == 0  # pi components
                strength = zeeman_strength(Jl, Ju, Ml, Mu)
                if strength > 0
                    push!(π_strength, strength)
                    push!(π_shift, gu*Mu - gl*Ml)
                end
            elseif ΔM == 1  # sigma^+ or b components
                strength = zeeman_strength(Jl, Ju, Ml, Mu)
                if strength > 0
                    push!(σb_strength, strength)
                    push!(σb_shift, gu*Mu - gl*Ml)
                end
            end
        end
    end
    return (
        SVector(σr_strength...), SVector(σr_shift...),
        SVector(π_strength...), SVector(π_shift...),
        SVector(σb_strength...), SVector(σb_shift...),
    )
end
