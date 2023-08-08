"""
Functions for line extinction
"""


"""
Doppler width for mass in kg, temperature in K
"""
function doppler_width(λ0, mass, temperature::T)::T where T <: AbstractFloat
    return λ0 / ustrip(c_0) * sqrt(2 * ustrip(k_B) * temperature / mass)
end


"""
Damping constant for γ in rad / s, λ and ΔλD in nm.
"""
function damping(γ::T, λ::Real, ΔλD::Real)::T where T <: Real
    c1 = 1 / ustrip((4 * π * c_0) |> u"nm / s")
    return c1 * γ * λ^2 / ΔλD
end


function calc_broadening(
    γ_params::LineBroadening,
    temperature::T,
    e_density::T,
    h_neutral_density::T
)::T where T <: AbstractFloat
    γ = zero(T)
    # Natural broadening
    γ += γ_params.natural
    # van der Waals broadening (ABO or Unsold recipes), perturber is H I
    γ += _γ_add(γ_params.hydrogen_const, γ_params.hydrogen_exp, temperature, h_neutral_density)
    # Quadratic Stark broadening, perturber are electrons
    γ += _γ_add(γ_params.electron_const, γ_params.electron_exp, temperature, e_density)
    # Linear Stark broadening
    if γ_params.stark_linear_const != 0.0
        γ += γ_params.stark_linear_const * e_density ^ γ_params.stark_linear_exp
    end
    return γ
end


function _γ_add(
    multiplier::AbstractVector{<: Real},
    exponent::AbstractVector{<: Real},
    base::AbstractFloat,
    perturber_density::T,
)::T where T <: AbstractFloat
    γ = zero(T)
    for (mult, expo) in zip(multiplier, exponent)
        γ += mult * base ^ expo
    end
    γ *= perturber_density
    return γ
end


function create_voigt_itp(a::AbstractArray{T}, v::AbstractRange) where {T <: AbstractFloat}
    tmp = Array{T, 2}(undef, length(a), length(v))
    Threads.@threads for iv in eachindex(v)
        for (ia, a_value) in enumerate(a)
            z = a_value - v[iv] * im
            tmp[ia, iv] = real(erfcx(z))
        end
    end
    return linear_interpolation((a, v), tmp, extrapolation_bc=Line())
end
