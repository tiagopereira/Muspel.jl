"""
Various utility functions.
"""

"""
Calculates the Blackbody (Planck) function per wavelength in nm
and temperature in K. Outputs in kW m^-2 nm^-1.
"""
function blackbody_λ(λ::T1, temperature::T2)::T2 where {T1 <: Real, T2 <: Real}
    mult = ustrip((2h * c_0^2) |>  u"kW * m^-2 * nm^4")
    return mult * λ^-5 / (exp(h * c_0 / k_B / (λ * temperature * u"K * nm")) - 1)
end


"""
    incline_data_inv!(
            data::AbstractArray{<: Real, 3},
            z::AbstractVector,
            dx::Real,
            dy::Real,
            μ::Real,
            ϕ::Real
    )

Transforms a 3D array into an inclined coordinate system, according to
a polar angle given by μ = cos(θ), and an azimuthal angle ϕ.
Uses cubic spline interpolation, and updates `data` in-place.
Based on trnslt.f90 from Åke Nordlund.

This version only works if height is the last dimension in `data`.
This is 3-5x faster than `incline_data!` because the loop order is optimised
for this memory configuration.
"""
function incline_data_inv!(
        data::AbstractArray{<: Real, 3},
        z::AbstractVector,
        dx::Real,
        dy::Real,
        μ::Real,
        ϕ::Real
)
    θ = acos(μ)
    sinθ = sin(θ)
    tanθ = sinθ / μ
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)
    ∂x∂z = tanθ * cosϕ
    ∂y∂z = tanθ * sinϕ
    nx, ny, nz = size(data)
    @assert nz == length(z)
    ε = 1.0e-6
    if abs(∂x∂z) > ε   # μ shift in the x dimension
        Threads.@threads for n in 1:nz
            shift_x = ∂x∂z * z[n] / (nx*dx)
            (k, ac, bc, ad, bd) = _spline_coeffs(shift_x, nx)
            buf = data[:, :, n]
            for l in 1:nx
                m1, p0, p1, p2 = _spline_stencil(l, k, nx)
                for m in 1:ny
                    data[l,m,n] = ac*buf[p0,m] + bc*buf[p1,m] - ad*buf[m1,m] + bd*buf[p2,m]
                end
            end
        end
    end
    if abs(∂y∂z) > ε   # ϕ shift in the y dimension
        Threads.@threads for n in 1:nz
            shift_y = ∂y∂z * z[n] / (ny*dy)
            (k, ac, bc, ad, bd) = _spline_coeffs(shift_y, ny)
            buf = data[:, :, n]
            for m in 1:ny
                m1, p0, p1, p2 = _spline_stencil(m, k, ny)
                for l in 1:nx
                    data[l,m,n] = ac*buf[l,p0] + bc*buf[l,p1] - ad*buf[l,m1] + bd*buf[l,p2]
                end
            end
        end
    end
end


"""
    incline_data!(
            data::AbstractArray{<: Real, 3},
            z::AbstractVector,
            dx::Real,
            dy::Real,
            μ::Real,
            ϕ::Real
    )

Transforms a 3D array into an inclined coordinate system, according to
a polar angle given by μ = cos(θ), and an azimuthal angle ϕ.
Uses cubic spline interpolation, and updates `data` in-place.

This version only works if height is the first dimension in `data`.
"""
function incline_data!(
    data::AbstractArray{<: Real, 3},
    z::AbstractVector,
    dx::Real,
    dy::Real,
    μ::Real,
    ϕ::Real
)
    data_inv = PermutedDimsArray(data, (3,2,1))
    incline_data_inv!(data_inv, z, dx, dy, μ, ϕ)
end


"""
Computes coefficients for cubic spline interpolation.
"""
function _spline_coeffs(shift, n)
    shift = mod(shift, 1)
    if shift < 0
        shift += 1
    end
    shift *= n
    k = floor(Int64, shift)
    p = shift - k
    k += n
    q = 1 - p
    af = q + p*q*(q - p)
    bf = p - p*q*(q - p)
    ad = p * q * q * 0.5f0
    bd = -p * q * p * 0.5f0
    ac = af - bd
    bc = bf + ad
    return (k, ac, bc, ad, bd)
end


"""
Gets stencil coordinates for cubic spline.
"""
function _spline_stencil(index, k, n)
    m1 = mod(index + k - 2, n) + 1
    p0 = mod(index + k - 1, n) + 1
    p1 = mod(index + k, n) + 1
    p2 = mod(index + k + 1, n) + 1
    return (m1, p0, p1, p2)
end
