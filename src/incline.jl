"""
Tools for inclining atmospheres.
"""

"""
    incline_atmos(atmos_in::AbstractAtmos3D, μ::Real, φ::Real)

Transforms a 3D atmosphere into an inclined coordinate system, given by a rotation
by a polar angle θ, given by μ = cos(θ) and an azimuthal angle φ. The output is
an atmosphere of the same type, where the quantities have been interpolated to the
same number of depth points and the vector quantities projected onto the new axes.

Assumes the following:

1. Order of the axes in 3D arrays is (z, y, x)
2. Horizontally periodic boundary conditions
3. The first index in the height direction is fixed in the polar rotation
4. A right handed system, so that:

```
      z ↑.            - θ is the clockwise polar rotation from the z axis
        |  .          - φ is the clockwise azimuthal rotation from the x axis
        |    .
        |      .
        |        *
        |      . .
        | θ  .   .
        |  .     .
        |.       .
        +--------.-----→ y
       /  .      .   /
      /  φ  .    .
     /         . . /
    /- - - - - - .
  x
```
"""
function incline_atmos(atmos_in::AbstractAtmos3D, μ::Real, φ::Real)
    atmos = deepcopy(atmos_in)
    dx = abs(atmos.x[2] - atmos.x[1])
    dy = abs(atmos.y[2] - atmos.y[1])
    for name in fieldnames(typeof(atmos))
        var = getfield(atmos, name)
        if ndims(var) == 3
           incline_data!(var, atmos.z, dx, dy, μ, φ)
        end
    end
    # adjust physical sizes
    atmos.z ./= μ
    atmos.x .*= μ * cos(φ)
    atmos.y .*= μ * sin(φ)
    # project vector quantities
    project_vector!(atmos.velocity_x, atmos.velocity_y, atmos.velocity_z, μ, φ)
    if all([f in fieldnames(typeof(atmos)) for f in [:Bx, :By, :Bz]])
        project_vector!(atmos.Bx, atmos.By, atmos.Bz, μ, φ)
    end
    return atmos
end


"""
    project_vector!(vx::A, vy::A, vz::A, μ::Real, φ::Real)

Projects a vector in 3D space according to the rotation by an polar angle
θ, given by μ = cos(θ) and an azimuthal angle φ. The inputs are the
vector components in the x, y, and z axes, for an array of a given size
(typically 3D).

Assumes the following:

1. θ is the clockwise polar rotation from the z axis
2. φ is the clockwise azimuthal rotation from the x axis
3. The rotation of the vector is the same as the rotation from the
   axes (x, y, z) into a new system (x', y', z') given by the rotation matrix:

```
               Polar          Azimuthal
   ⎡x'⎤ = ⎡cosθ  0  -sinθ⎤⎡ cosφ  sinφ  0⎤⎡x⎤
   ⎜y'⎥   ⎜  0   1    0  ⎥⎜-sinφ  cosφ  0⎥⎜y⎥
   ⎣z'⎦   ⎣sinθ  0   cosθ⎦⎣   0     0   1⎦⎣z⎦

   ⌈x'⎤ = ⎡cosθcosφ  cosθsinφ  -sinφ⎤⎡x⎤
   |y'⎥   ⎜ -sinφ      cosφ       0 ⎥⎜y⎥
   ⌊z'⎦   ⎣sinθcosφ  sinθsinφ   cosθ⎦⎣z⎦
```
"""
function project_vector!(vx::A, vy::A, vz::A, μ::Real, φ::Real) where A <: AbstractArray{<:Real}
    cosθ = μ
    sinθ = sqrt(1 - μ^2)
    cosφ = cos(φ)
    sinφ = sin(φ)
    Threads.@threads for i in eachindex(vx)
        proj_x = vx[i]*cosθ*cosφ + vy[i]*cosθ*sinφ - vz[i]*sinθ
        proj_y = -vx[i]*sinφ + vy[i]*cosφ
        proj_z = vx[i]*sinθ*cosφ + vy[i]*sinθ*sinφ + vz[i]*cosθ
        vx[i] = proj_x
        vy[i] = proj_y
        vz[i] = proj_z
    end
    return nothing
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
    # Invert θ so that inclination goes towards increasing x or y values
    μ = -μ
    θ = acos(μ)
    sinθ = sin(θ)
    tanθ = sinθ / μ
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)
    ∂x∂z = tanθ * cosϕ
    ∂y∂z = tanθ * sinϕ
    nx, ny, nz = size(data)
    @assert nz == length(z)
    @assert dx != 0
    @assert dy != 0
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
