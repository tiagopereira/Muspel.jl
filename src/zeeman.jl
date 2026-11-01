"""
Calculates the effective Landé factor.
"""
function g_eff(g1, g2, J1, J2)
    d = J1 * (J1 + 1) - J2 * (J2 + 1)
    return (g1 + g2) * 0.5f0 + (g1 - g2) * d * 0.25f0
end


"""
Utility function to compute Landé factors from various types of coupling.
From pp. 76-77 of
[Landi Degl'Innocenti & Landolfi (2004)](https://ui.adsabs.harvard.edu/abs/2004ASSL..307.....L/abstract).
"""
function γ_g(A, B, C)
    if A == 0
        return 0f0
    else
        return 0.5f0 * (A * (A + 1) + B * (B + 1) - C * (C + 1)) / (A * (A + 1))
    end
end


"""
Calculates Landé factor from LS coupling.
"""
g_LS(S, L, J) = 1f0 + γ_g(J, S, L)


"""
Calculates Landé factor from jj coupling.
"""
g_jj(J, j1, j2, l1, l2) = 1f0 + γ_g(J, j1, j2)*γ_g(j1, 1//2, l1) + γ_g(J, j1, j2)*γ_g(j2, 1//2, l2)


"""
Calculates Landé factor from JK coupling.
"""
g_JK(J, K, J1, S1, L1, l) = 1f0 + γ_g(J, 1//2, K) + γ_g(J, K, 1//2)*γ_g(K, J1, l)*γ_g(J1, S1, L1)


"""
Calculates strengths of Zeeman components. These strengths are normalised.
"""
function zeeman_strength(J_l, J_u, M_l, M_u)
    J, M = J_l, M_l
    if J_l == J_u == 0
        throw(DomainError("Invalid transition, J_u = J_l = 0"))
    end
    if M_u - M_l == 1
        if J_u - J_l == 1
            strength = (3 * (J + M + 1) * (J + M + 2)) / (2 * (J + 1) * (2 * J + 1) * (2 * J + 3))
        elseif J_u - J_l == 0
            strength = (3 * (J - M) * (J + M + 1)) / (2 * J * (J + 1) * (2 * J + 1))
        elseif J_u - J_l == -1
            strength = (3 * (J - M) * (J - M - 1)) / (2 * J * (2 * J - 1) * (2 * J + 1))
        else
            throw(DomainError("Invalid transition, J_u - J_l != -1, 0, 1"))
        end
    elseif M_u - M_l == 0
        if J_u - J_l == 1
            strength = (3 * (J - M + 1) * (J + M + 1)) / ((J + 1) * (2 * J + 1) * (2 * J + 3))
        elseif J_u - J_l == 0
            strength = 3 * M^2 / (J * (J + 1) * (2 * J + 1))
        elseif J_u - J_l == -1
            strength = (3 * (J - M) * (J + M)) / (J * (2 * J - 1) * (2 * J + 1))
        else
            throw(DomainError("Invalid transition, J_u - J_l != -1, 0, 1"))
        end
    elseif M_u - M_l == -1
        if J_u - J_l == 1
            strength = (3 * (J - M + 1) * (J - M + 2)) / (2 * (J + 1) * (2 * J + 1) * (2 * J + 3))
        elseif J_u - J_l == 0
            strength = (3 * (J + M) * (J - M + 1)) / (2 * J * (J + 1) * (2 * J + 1))
        elseif J_u - J_l == -1
            strength = (3 * (J + M) * (J + M - 1)) / (2 * J * (2 * J - 1) * (2 * J + 1))
        else
            throw(DomainError("Invalid transition, J_u - J_l != -1, 0, 1"))
        end
    else
        throw(DomainError("Invalid transition, M_u - M_l != -1, 0, 1"))
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
    if count("E", label) > 0 && count("O", label) > 0
        if findlast("E", label) > findlast("O", label)
            label = split(label, "E")[end-1]
        else
            label = split(label, "O")[end-1]
        end
    elseif count("E", label) > 0
        label = split(label, "E")[end-1]
    elseif count("O", label) > 0
        label = split(label, "O")[end-1]
    end
    tmp = split(label)
    if length(tmp) == 0
        throw(ArgumentError("Found empty label"))
    else
        term = tmp[end]
        try
            multiplicity = parse(Int, term[end-1])
            S = (multiplicity - 1) // 2
            L = orbital[uppercase(term[end])]
            return (S, L)
        catch e
            throw(ArgumentError("Could not parse term symbol: $term"))
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


"""
    field_geometry(Bx, By, Bz)

Computes the geometry of a magnetic field vector needed for the Zeeman profiles,
without any trigonometric functions. `z` is the line of sight, so `γ` is the
inclination and `χ` the azimuth (from `x`) of the field.

Returns `(B, cosγ, sin²γ, cos2χ, sin2χ)`. For `B = 0` or a field along the line of sight
(where `χ` is undefined) it returns `sin²γ = 0` and a fixed, arbitrary azimuth.
"""
@inline function field_geometry(Bx::T, By::T, Bz::T) where T <: AbstractFloat
    B⊥² = Bx^2 + By^2
    B = sqrt(B⊥² + Bz^2)
    if B == 0
        return zero(T), one(T), zero(T), one(T), zero(T)
    end
    cosγ = Bz / B
    sin²γ = B⊥² / B^2
    if B⊥² == 0
        cos2χ, sin2χ = one(T), zero(T)
    else
        cos2χ = (Bx^2 - By^2) / B⊥²
        sin2χ = 2 * Bx * By / B⊥²
    end
    return B, cosγ, sin²γ, cos2χ, sin2χ
end


"""
    stokes_profiles(line, voigt_itp, a, v, ΔλD, ΔλB, cosγ, sin²γ, cos2χ, sin2χ)

Computes line profiles under Zeeman effect, for all
Stokes parameters and also including magneto-optical effects.
`ΔλB` is the Zeeman wavelength shift, and the field geometry is given by
`cosγ`, `sin²γ`, `cos2χ` and `sin2χ` (see [`field_geometry`](@ref)).

Returns a 7-element tuple with: (ϕI, ϕQ, ϕU, ϕV, ψQ, ψU, ψV).

Takes into account the symmetry of H and anti-symmetry of L,
so voigt_itp only needs to cover positive values of v.
"""
@inline function stokes_profiles(
    line::AtomicLine,
    voigt_itp::Interpolations.AbstractInterpolation{<:Number, 2},
    a::T,
    v::T,
    ΔλD::T,
    ΔλB::T,
    cosγ::T,
    sin²γ::T,
    cos2χ::T,
    sin2χ::T,
)::NTuple{7, T} where T <: AbstractFloat
    ΔλB_D = ΔλB / ΔλD
    ϕ0, ψ0 = _zeeman_group(voigt_itp, a, v, ΔλB_D, line.π_strength, line.π_shift)  # ΔMj = 0
    ϕm, ψm = _zeeman_group(voigt_itp, a, v, ΔλB_D, line.σr_strength, line.σr_shift)  # ΔMj = -1
    ϕp, ψp = _zeeman_group(voigt_itp, a, v, ΔλB_D, line.σb_strength, line.σb_shift)  # ΔMj = +1
    norm = one(T) / (sqrt(T(π)) * ΔλD)
    ϕ0, ϕm, ϕp = ϕ0 * norm, ϕm * norm, ϕp * norm
    ψ0, ψm, ψp = ψ0 * norm, ψm * norm, ψp * norm
    ϕΔ = (ϕ0 - (ϕp + ϕm) / 2) / 2
    ϕI = ϕΔ * sin²γ + (ϕp + ϕm) / 2
    ϕQ = ϕΔ * sin²γ * cos2χ
    ϕU = ϕΔ * sin²γ * sin2χ
    ϕV = (ϕm - ϕp) * cosγ / 2
    ψΔ = (ψ0 - (ψp + ψm) / 2) / 2
    ψQ = ψΔ * sin²γ * cos2χ
    ψU = ψΔ * sin²γ * sin2χ
    ψV = (ψm - ψp) * cosγ / 2
    return (ϕI, ϕQ, ϕU, ϕV, ψQ, ψU, ψV)
end


"""
Sums the strength-weighted Voigt (real) and Faraday-Voigt (imaginary) profiles of a group
of Zeeman components, without the `1 / (√π ΔλD)` normalisation. `ΔλB_D` is `ΔλB / ΔλD`.
"""
@inline function _zeeman_group(voigt_itp, a::T, v::T, ΔλB_D::T, strength, shift) where T
    ϕ = ψ = zero(T)
    for i in eachindex(strength)
        v_zeeman = v + shift[i] * ΔλB_D
        voigt = voigt_itp(a, abs(v_zeeman))
        ϕ += strength[i] * real(voigt)
        ψ += copysign(strength[i], v_zeeman) * imag(voigt)
    end
    return ϕ, ψ
end


"""
Builds the reduced 4x4 Stokes propagation matrix, defined as:

```
A = α_l / α_I ⎡0    ϕQ   ϕU   ϕV⎤
              ⎜ϕQ    0   ψV  -ψU⎥
              ⎜ϕU  -ψV    0   ψQ⎥
              ⎣ϕV   ψU  -ψQ    0⎦

```
where `α_I = α_cont + α_line * ϕI`.
"""
function StokesA(ϕQ, ϕU, ϕV, ψQ, ψU, ψV, α_line, α_I)
    α = α_line / α_I
    return @SMatrix[    0  α*ϕQ  α*ϕU  α*ϕV;
                     α*ϕQ     0  α*ψV -α*ψU;
                     α*ϕU -α*ψV     0  α*ψQ;
                     α*ϕV  α*ψU -α*ψQ     0]
end


"""
Coefficient `e / (4π m_e c)` in nm⁻¹ T⁻¹, such that the Larmor wavelength shift
in nm is `ΔλB = larmor_coeff * λ0² * B`, with `λ0` in nm and `B` in T.
"""
const larmor_coeff = ustrip((e / (4 * π * m_e * c_0)) |> u"nm^-1 * T^-1")
