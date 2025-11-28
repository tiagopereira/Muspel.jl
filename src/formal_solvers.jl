"""
Tools for computing the formal solution of the radiative transfer equation.
"""


"""
    piecewise_1D_linear(
        z::AbstractVector{T},
        α::AbstractVector{T},
        source_function::AbstractVector{T};
        to_end::Bool=false,
        initial_condition=:source
    ) where T <: AbstractFloat

DEPRECATED, use piecewise_1D_linear!
Compute piecewise integration of the radiative transfer equation,
assuming linear integration of the source function, for a given
height `z`, extinction `α` and `source_function`. The optional
keyword argument `to_end` defines the direction of the integration:
if `false` (default) will start to integrate intensity from the last
element to the first, and if `true` will integrate from the first
element to the last. `initial_condition` can take two values: `:zero` for
no radiation, or `:source` (default) to take the source function at the
starting point.
"""
function piecewise_1D_linear(
    z::AbstractVector{T},
    α::AbstractVector{T},
    source_function::AbstractVector{T};
    to_end::Bool=false,
    initial_condition=:source
) where T <: AbstractFloat
    ndep = length(z)
    if to_end
        start = 1
        incr = 1
        depth_range = 2:ndep
    else
        start = ndep
        incr = -1
        depth_range = ndep-1:-1:1
    end
    intensity = similar(source_function)
    if initial_condition == :source
        intensity[start] = source_function[start]
    elseif initial_condition == :zero
        intensity[start] *= 0
    else
        throw(ErrorException("NotImplemented initial condition $initial_condition"))
    end
    for i in depth_range
        Δτ = abs(z[i] - z[i-incr]) * (α[i] + α[i-incr]) / 2
        w1, w2, w3 = _w3(Δτ)
        intensity[i] = w1*intensity[i-incr] + w2*source_function[i-incr] + w3*source_function[i]
    end
    return intensity
end


"""
    piecewise_1D_linear!(
        z::AbstractVector{T},
        α::AbstractVector{T},
        source_function::AbstractVector{T},
        intensity::AbstractVector{T};
        to_end::Bool=false,
        initial_condition=:source
    ) where T <: AbstractFloat

Compute piecewise integration of the radiative transfer equation,
assuming linear integration of the source function, for a given
height `z`, extinction `α` and `source_function` and `intensity`
(existing array where the output will be saved into). The optional
keyword argument `to_end` defines the direction of the integration:
if `false` (default) will start to integrate intensity from the last
element to the first, and if `true` will integrate from the first
element to the last. `initial_condition` can take two values: `:zero` for
no radiation, or `:source` (default) to take the source function at the
starting point.
"""
function piecewise_1D_linear!(
    z::AbstractVector{T},
    α::AbstractVector{T},
    source_function::AbstractVector{T},
    intensity::AbstractVector{T};
    to_end::Bool=false,
    initial_condition=:source
) where T <: AbstractFloat
    ndep = length(z)
    if to_end
        start = 1
        incr = 1
        depth_range = 2:ndep
    else
        start = ndep
        incr = -1
        depth_range = ndep-1:-1:1
    end
    if initial_condition == :source
        intensity[start] = source_function[start]
    elseif initial_condition == :zero
        intensity[start] *= 0
    else
        throw(ErrorException("NotImplemented initial condition $initial_condition"))
    end
    for i in depth_range
        Δτ = abs(z[i] - z[i-incr]) * (α[i] + α[i-incr]) / 2
        w1, w2, w3 = _w3(Δτ)
        intensity[i] = w1*intensity[i-incr] + w2*source_function[i-incr] + w3*source_function[i]
    end
    return nothing
end


"""
    function feautrier(
        z::Array{<:Unitful.Length{T}, 1},
        α::Array{<:PerLength{T}, 1},
        source_function::Array{<:Unitful.Quantity, 1}
    ) where T <: AbstractFloat

Calculate solution to radiative transfer equation using the Feautrier method.
Uses algorithm from [Rybicki & Hummer, 1991, A&A 245](https://ui.adsabs.harvard.edu/abs/1991A&A...245..171R).
Returns the height-dependent Feautrier variable `P`:

```math
P \\equiv 1/2 (I^+ + I^-)
```

Currently operates under the following assumptions:
* The first index of the variables is the top of the atmosphere
* The boundary conditions are zero radiation at the top and source function at the bottom

Therefore, the emergent intensity is `2 * P[1]`, since ``I^-[1]=0``.

Not properly tested, use with care!
"""
function feautrier(
    z::AbstractVector{T},
    α::AbstractVector{T},
    source_function::AbstractVector{T},
) where T <: AbstractFloat
    ndep = length(z)
    F = Array{T}(undef, ndep)
    Z = similar(source_function)
    P = similar(source_function)
    Δτ = Array{T}(undef, ndep)
    Δτ[end] = zero(T)
    for i in 1:ndep-1
        Δτ[i] = (z[i] - z[i+1]) * (α[i] + α[i+1]) / 2
    end
    H1 = one(T) + 2 / Δτ[1]
    C1 = 2 / Δτ[1]^2
    Hn = one(T) + 2 / Δτ[end - 1]
    An = 2 / Δτ[end - 1]^2
    # Start elimination
    F[1] = H1 / C1
    Z[1] = source_function[1] / (H1 + C1)
    for i in 2:ndep-1
        Δτ_mid = (Δτ[i] + Δτ[i - 1]) / 2
        A = 1 / (Δτ_mid * Δτ[i - 1])
        C = 1 / (Δτ_mid * Δτ[i])
        F[i] = (one(T) + A * F[i - 1] / (one(T) + F[i - 1])) / C
        Z[i] = (source_function[i] + A * Z[i - 1]) / (C * (one(T) + F[i]))
    end
    # Now backsubstitution
    P[end] = (source_function[end] + An * Z[end - 1]) /
        (Hn + An * (F[end - 1] / (one(T) + F[end - 1])))
    for i in ndep-1:-1:1
        P[i] = P[i + 1] / (one(T) + F[i]) + Z[i]
    end
    return P
end


"""
Computes weights for linear integration of source function,
approximating `exp(-Δτ)` for very small and very large values of `Δτ`.
From Jaime de la Cruz Rodriguez, adapted from RH.
"""
function _w3(Δτ::T) where T <: AbstractFloat
    if Δτ > 40   # assume exp(-Δτ) = 0
        w1 = zero(T)
        w2 = one(T) / Δτ
        w3 = one(T) - w2
    elseif Δτ > 0.01   # normal range
        w1 = exp(-Δτ)
        u0 = (one(T) - w1) / Δτ
        w2 = u0 - w1
        w3 = one(T) - u0
    else   # Taylor expansion at Δτ=0
        w1 = one(T) - Δτ + Δτ^2 / 2
        w2 = (one(T)/2 - Δτ/3) * Δτ
        w3 = (one(T)/2 - Δτ/6) * Δτ
    end
    return w1, w2, w3
end


#=
Piecewise cubic Bezier based on de la Cruz Rodriguez & Piskunov (2013),
https://ui.adsabs.harvard.edu/abs/2013ApJ...764...33D

Adapted from RH / STiC.
=#

"""
Compute interpolation coefficients for cubic Bézier interpolation
for a given optical depth interval δ and following
[de la Cruz Rodriguez & Piskunov (2013)](https://ui.adsabs.harvard.edu/abs/2013ApJ...764...33D),
after eq. (20).

Uses approximations for low or high optical depths to avoid
numerical instabilities.
"""
function bezier3_coeffs(δin::T)::NTuple{5, T} where T <: AbstractFloat
    δ = Float64(δin)  # Becomes imprecise in Float32, need to ensure Float64
    if δ < 0.05
        ϵ = 1. - δ + δ^2 / 2 - δ^3 / 6
        α = δ / 4 - δ^2 / 5 + δ^3 / 12
        β = δ / 4 - δ^2 / 20 + δ^3 / 120
        γ = δ / 4 - δ^2 * 3/20 + δ^3 / 20
        θ = δ / 4 - δ^2 / 10 + δ^3 / 40
    elseif δ > 30
        ϵ = 0.
        α = 6 / δ^3
        β = 1. - α*(1-δ) - 3/δ
        γ = α * (δ-3)
        θ = 3 * (α + (δ-4)/δ^2)
    else
        ϵ = exp(-δ)
        α = -(-6 + (6 + 6*δ + 3*δ^2 + δ^3)*ϵ) / δ^3
        β = (-6 + δ*(6 + δ * (δ-3)) + 6*ϵ) / δ^3
        γ = 3 * (2*δ - 6 + ϵ*(6 + δ * (δ+4))) / δ^3
        θ = 3 * (-2 * ϵ * (δ+3) + 6 + (δ-4)*δ) / δ^3
    end
    return (convert(T, α), convert(T, β), convert(T, γ), convert(T, θ), convert(T, ϵ))
end


"""
Compute derivative of a monotonic interpolation function following
[Steffen (1990)](https://ui.adsabs.harvard.edu/abs/1990A%26A...239..443S),
eq. (11).
"""
function cent_deriv(odx, dx, yu, y0, yd)
    S0 = (yd - y0) / dx
    Su = (y0 - yu) / odx
    P0 = abs((Su*dx + S0*odx) / (odx+dx)) / 2
    return (sign(S0) + sign(Su)) * min(abs(Su), min(abs(S0), P0))
end


"""
    function piecewise_1D_bezier3!(
        z::AbstractVector{T},
        α::AbstractVector{T},
        source::AbstractVector{T},  # source function
        intensity::AbstractVector{T};
        to_end::Bool=false,
        initial_condition=:source
    ) where T <: AbstractFloat

Compute piecewise integration of the radiative transfer equation,
assuming a cubic Bézier integration of the source function, following
[de la Cruz Rodriguez & Piskunov (2013)](https://ui.adsabs.harvard.edu/abs/2013ApJ...764...33D),
and adapted from the implementation in RH / STiC.
Calculates for a given height `z`, extinction `α` and `source_function`
and `intensity` (existing array where the output will be saved into).
The optional keyword argument `to_end` defines the direction of the
integration: if `false` (default) will start to integrate intensity from
the last element to the first, and if `true` will integrate from the first
element to the last. `initial_condition` can take two values: `:zero` for
no radiation, or `:source` (default) to take the source function at the
starting point.
"""
function piecewise_1D_bezier3!(
    z::AbstractVector{T},
    α::AbstractVector{T},
    source::AbstractVector{T},  # source function
    intensity::AbstractVector{T};
    to_end::Bool=false,
    initial_condition=:source
) where T <: AbstractFloat
    ndep = length(z)
    if to_end
        start = 1
        incr = 1
        depth_range = 2:ndep
        iend = ndep
    else  # to_obs in RH, default
        start = ndep
        incr = -1
        depth_range = ndep-1:-1:1
        iend = 1
    end
    if initial_condition == :source
        intensity[start] = source[start]
    elseif initial_condition == :zero
        intensity[start] = zero(T)
    else
        throw(ErrorException("NotImplemented initial condition $initial_condition"))
    end

    # set variables for first iteration to allow simple
    # shift for all next iterations
    i = start + incr
    dsup = abs(z[i] - z[i-incr])
    dsdn = abs(z[i+incr] - z[i])
    dα_up = (α[i] - α[i-incr]) / dsup
    dα_dn = zero(T)
    dτ_dw = zero(T)

    #  dα/ds at central point
    dα_c = cent_deriv(dsup, dsdn, α[i-incr], α[i], α[i+incr])
    # upwind path_length (Bezier3 integration)
    c1 = α[i] - dsup/3 * dα_c
    c2 = α[i-incr] + dsup/3 * dα_up
    dτ_uw = dsup * (α[i] + α[i-incr] + c1 + c2) / 4

    # dS/dtau at upwind point
    dS_up = (source[i] - source[i-incr]) / dτ_uw

    for i in depth_range
        if i != iend
            # downwind path length
            dsdn = abs(z[i] - z[i-incr])
            # dα/ds at downwind point
            if abs(i - iend) > 1
                dsdn2 = abs(z[i + 2*incr] - z[i+incr])
                dα_dn = cent_deriv(dsdn, dsdn2, α[i], α[i+incr], α[i + 2*incr])
            else
                dα_dn = (α[i+incr] - α[i]) / dsdn
            end
            #  Do not clip control points, it fails if there is much stimulated emission
            c1 = α[i] + dsdn/3 * dα_c
            c2 = α[i+incr] - dsdn/3 * dα_dn

            # downwind optical path length
            dτ_dw = dsdn * (α[i] + α[i+incr] + c1 + c2) / 4

            #  dS/dt at central point
            dS_c = cent_deriv(dτ_uw, dτ_dw, source[i-incr], source[i], source[i+incr])
        else
            # Last interval, use linear approx for the derivative at
	        # the central point
             dS_c = (source[i] - source[i-incr]) / dτ_uw
        end

        # compute interpolation parameters
        (χ, β, γ, θ, ϵ) = bezier3_coeffs(dτ_uw)

        # Source function control points
        c1 = source[i] - dS_c * dτ_uw / 3
        c2 = source[i-incr] + dS_up * dτ_uw / 3

        #  Solve integral in this interval
        intensity[i] = ϵ*intensity[i-incr] + β*source[i] + χ*source[i-incr] + θ*c1 + γ*c2

        # Re-use downwind quantities for next upwind position
        dsup = dsdn
        dα_up = dα_c
        dα_c = dα_dn
        dτ_uw = dτ_dw
        dS_up = dS_c
    end
    return nothing
end
