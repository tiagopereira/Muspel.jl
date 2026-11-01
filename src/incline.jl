"""
Tools for inclining atmospheres.
"""

"""
    incline_atmos(
        atmos_in::AbstractAtmosphere{3},
        μ::Real,
        φ::Real;
        interpolation=:steffen,
        vector_interpolation=:weno4,
    )

Transforms a 3D atmosphere into an inclined coordinate system, given by a rotation
by a polar angle θ, given by μ = cos(θ) and an azimuthal angle φ. The output is
an atmosphere of the same type, where the quantities have been interpolated to the
same number of depth points and the vector quantities projected onto the new axes.

Scalar quantities (temperature and densities, which cannot be negative) are
interpolated with `interpolation`, and vector quantities (velocity and magnetic field)
with `vector_interpolation`. See [`incline_data!`](@ref) for the available methods.
The default `:steffen` is monotonic, so it does not create new extrema and keeps
positive quantities positive, while `:weno4` is more accurate in smooth regions but
only essentially non-oscillatory.

Assumes the following:

1. Order of the axes in 3D arrays is (z, y, x)
2. Horizontally periodic boundary conditions
3. The first index in the height direction is fixed in the polar rotation
4. A right handed system, so that:

```
      z  ↑.           - θ is the clockwise polar rotation from the z axis
         |  .         - φ is the clockwise azimuthal rotation from the x axis
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
     /_ _ _ _ _ _ .
    /
  x
```
"""
function incline_atmos(
    atmos_in::AbstractAtmosphere{3},
    μ::Real,
    φ::Real;
    interpolation=:steffen,
    vector_interpolation=:weno4,
)
    atmos = deepcopy(atmos_in)
    incline_atmos!(atmos, atmos_in, μ, φ; interpolation, vector_interpolation)
    return atmos
end

"""
    incline_atmos!(
        atmos::AbstractAtmosphere{3},
        atmos_in::AbstractAtmosphere{3},
        μ::Real,
        φ::Real;
        interpolation=:steffen,
        vector_interpolation=:weno4,
        buffer=nothing,
    )

In-place version of [`incline_atmos`](@ref): the inclined atmosphere is written into
`atmos`, which must have the same type and size as `atmos_in` and not share any arrays
with it. Repeated calls (e.g. for many directions) avoid allocating a new atmosphere.
`buffer` can be a preallocated array with the same size as the data, used as workspace
when the inclination is along both horizontal axes.
"""
function incline_atmos!(
    atmos::AbstractAtmosphere{3},
    atmos_in::AbstractAtmosphere{3},
    μ::Real,
    φ::Real;
    interpolation=:steffen,
    vector_interpolation=:weno4,
    buffer=nothing,
)
    if keys(atmos_in.velocity) != (:x, :y, :z)
        throw(ArgumentError(
            "incline_atmos requires an atmosphere with x, y, z velocity components"
        ))
    end
    if typeof(atmos) != typeof(atmos_in) || size(atmos) != size(atmos_in)
        throw(ArgumentError("atmos and atmos_in must have the same type and size"))
    end
    if atmos.temperature === atmos_in.temperature
        throw(ArgumentError("atmos and atmos_in must not share arrays"))
    end
    if buffer === nothing
        buffer = similar(atmos_in.temperature)
    end
    dx = abs(atmos_in.x[2] - atmos_in.x[1])
    dy = abs(atmos_in.y[2] - atmos_in.y[1])
    nz, ny, nx = size(atmos_in)
    T = eltype(atmos_in.temperature)
    scalar_scheme = _interp_scheme(interpolation)
    vector_scheme = _interp_scheme(vector_interpolation)
    # The shifts are the same for all quantities that use the same interpolation
    scalar_plan = _incline_plan(scalar_scheme, T, nz, ny, nx, atmos_in.z, dx, dy, μ, φ)
    vector_plan = _incline_plan(vector_scheme, T, nz, ny, nx, atmos_in.z, dx, dy, μ, φ)
    for name in fieldnames(typeof(atmos_in))
        var_in = getfield(atmos_in, name)
        var = getfield(atmos, name)
        if var_in isa NamedTuple  # vector quantities
            for k in keys(var_in)
                _incline!(var_in[k], var[k], vector_scheme, vector_plan, buffer)
            end
        elseif var_in isa AbstractArray{<:Real, 3}  # scalar quantities
            _incline!(var_in, var, scalar_scheme, scalar_plan, buffer)
        end
    end
    # adjust physical sizes
    atmos.z .= atmos_in.z ./ μ
    atmos.x .= atmos_in.x .* (μ * cos(φ))
    atmos.y .= atmos_in.y .* (μ * sin(φ))
    # project vector quantities
    project_vector!(atmos.velocity.x, atmos.velocity.y, atmos.velocity.z, μ, φ)
    if has_magnetic_field(atmos)
        project_vector!(atmos.magnetic_field.x, atmos.magnetic_field.y,
                        atmos.magnetic_field.z, μ, φ)
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
3. The rotation of the vector is the same as the (inverse) rotation from the
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
    incline_data!(
            data_in::AbstractArray{<: Real, 3},
            data_out::AbstractArray{<: Real, 3},
            z::AbstractVector,
            dx::Real,
            dy::Real,
            μ::Real,
            ϕ::Real;
            interpolation=:steffen,
            buffer=nothing,
    )

Transforms a 3D array into an inclined coordinate system, according to
a polar angle given by μ = cos(θ), and an azimuthal angle ϕ.
Writes output into `data_out`, using the interpolation method given by `interpolation`:

* `:linear`: piecewise linear (second order, monotonic).
* `:cubic`: cubic Hermite spline with centred derivatives (third order). It can overshoot
  next to sharp gradients, and make positive quantities negative.
* `:steffen`: cubic Hermite with limited derivatives (Steffen 1990). It is monotonic
  between every pair of grid points, so it never overshoots, and the result is always
  bounded by the neighbouring values. Third order, except at local extrema.
* `:weno4`: fourth-order weighted essentially non-oscillatory (WENO) interpolation, using
  the smoothness indicators of [Janett et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019A%26A...624A.104J/abstract)
  (see also `WENO4.jl`), specialised for a uniform periodic grid. The smoothness
  indicators are scaled by the local magnitude of the data, so the result does not depend
  on units. It is the most accurate for smooth data and hardly overshoots at jumps
  (to ~1e-7 of the jump), but this is not strictly guaranteed.

The inclination is done as a sequence of periodic shifts of each horizontal plane
along x and then along y. A preallocated `buffer` with the same size as `data_in`
can be passed to avoid allocating a work array when shifting along both axes.
`data_in` and `data_out` must not be the same array.

Height must be the first dimension in `data_in` and `data_out`.
"""
function incline_data!(
        data_in::AbstractArray{<: Real, 3},
        data_out::AbstractArray{<: Real, 3},
        z::AbstractVector,
        dx::Real,
        dy::Real,
        μ::Real,
        ϕ::Real;
        interpolation=:steffen,
        buffer=nothing,
)
    scheme = _interp_scheme(interpolation)
    nz, ny, nx = size(data_in)
    plan = _incline_plan(scheme, eltype(data_in), nz, ny, nx, z, dx, dy, μ, ϕ)
    _incline!(data_in, data_out, scheme, plan, buffer)
    return nothing
end

# Interpolation of periodically shifted planes.
#
# Every z plane is shifted along x and/or y by a distance that is constant
# in the plane, so that the interpolation is that of a periodic uniform grid with a
# constant fractional offset `t` (between 0 and 1) and integer offset `k`.
# The output at index `l` is the interpolation between the input at `l + k` and `l + k + 1`,
# using a stencil of four input points (`ym`, `y0`, `y1`, `y2`).
# The `_interp` methods are pure scalar functions of the stencil values and of the
# per-plane coefficients, so they can be used unchanged inside a GPU kernel.

struct LinearInterp end
struct CubicInterp end
struct SteffenInterp end
struct Weno4Interp end

function _interp_scheme(name::Symbol)
    name === :linear && return LinearInterp()
    name === :cubic && return CubicInterp()
    name === :steffen && return SteffenInterp()
    name === :weno4 && return Weno4Interp()
    throw(ArgumentError(
        "Interpolation method :$name not recognised. Use :linear, :cubic, :steffen or :weno4."
    ))
end

# ε regularises the nonlinear weights of WENO4, relative to the (scaled) data
const WENO4_ε = 1.0f-6

# Size of the blocks of heights processed together when shifting
const SHIFT_BLOCK = 32

# Padding between the work arrays of different threads, in elements, so that threads do not
# share cache lines (false sharing slows down the shifting by around 50%)
const SHIFT_PAD = 32

"""
Shifts needed to incline data of size `(nz, ny, nx)` and type `T`, for the polar angle
given by μ = cos(θ) and the azimuthal angle ϕ: for each axis the integer shift of each
height `k` and the coefficients of the interpolation `scheme`, and whether any shift is
needed. It only depends on the geometry, so it can be reused for all quantities.
"""
function _incline_plan(scheme, ::Type{T}, nz, ny, nx, z, dx, dy, μ, ϕ) where T
    # Invert θ so that inclination goes towards increasing x or y values
    μ = -μ
    @assert abs(μ) > 0 "μ must be non-zero"
    θ = acos(μ)
    sinθ = sin(θ)
    tanθ = sinθ / μ
    cosϕ = cos(ϕ)
    sinϕ = sin(ϕ)
    ∂x∂z = tanθ * cosϕ
    ∂y∂z = tanθ * sinϕ
    @assert nz == length(z)
    @assert dx != 0
    @assert dy != 0
    ε = 1.0e-6
    kx, cx = _plane_shifts(scheme, ∂x∂z, z, nx*dx, nx, T)
    ky, cy = _plane_shifts(scheme, ∂y∂z, z, ny*dy, ny, T)
    return (shift_x=abs(∂x∂z) > ε, shift_y=abs(∂y∂z) > ε, kx=kx, cx=cx, ky=ky, cy=cy)
end

"""
Inclines `data_in` into `data_out` following `plan` (from `_incline_plan`). Shifting along
both axes needs a work array, which is `buffer` if given.
"""
function _incline!(data_in, data_out, scheme, plan, buffer)
    @assert size(data_out) == size(data_in)
    if !plan.shift_x && !plan.shift_y
        data_out .= data_in
    elseif plan.shift_x && plan.shift_y
        if buffer === nothing
            buffer = similar(data_in)
        end
        @assert size(buffer) == size(data_in)
        _shift_x!(buffer, data_in, scheme, plan)
        _shift_y!(data_out, buffer, scheme, plan)
    elseif plan.shift_x
        _shift_x!(data_out, data_in, scheme, plan)
    else
        _shift_y!(data_out, data_in, scheme, plan)
    end
    return nothing
end

_shift_x!(data_out, data_in, scheme, plan) =
    _shift_axis!(data_out, data_in, scheme, plan.kx, plan.cx)

# Same as shifting along x, with the horizontal axes swapped
_shift_y!(data_out, data_in, scheme, plan) =
    _shift_axis!(PermutedDimsArray(data_out, (1, 3, 2)),
                 PermutedDimsArray(data_in, (1, 3, 2)), scheme, plan.ky, plan.cy)

"""
Per-plane integer shift `k` (between 0 and `n - 1`) and the coefficients used by
the interpolation scheme, for a periodic axis with `n` points and total length `L`.
`∂∂z` is the horizontal displacement per unit of height, in the units of `L`.
"""
function _plane_shifts(scheme, ∂∂z, z, L, n, ::Type{T}) where T
    nz = length(z)
    k = Vector{Int}(undef, nz)
    t = Vector{T}(undef, nz)
    for i in 1:nz
        shift = mod(∂∂z * z[i] / L, 1) * n
        ki = floor(Int, shift)
        t[i] = T(shift - ki)
        k[i] = mod(ki, n)
    end
    return k, _coeffs(scheme, t)
end

_coeffs(::LinearInterp, t) = (t,)
_coeffs(::SteffenInterp, t) = (t,)
_coeffs(::Weno4Interp, t::Vector{T}) where T = (t, [x * (x - 1) / 2 for x in t])
function _coeffs(::CubicInterp, t::Vector{T}) where T
    ac = similar(t)
    bc = similar(t)
    ad = similar(t)
    bd = similar(t)
    for i in eachindex(t)
        p = t[i]
        q = 1 - p
        af = q + p*q*(q - p)
        bf = p - p*q*(q - p)
        ad[i] = p * q * q * T(0.5)
        bd[i] = -p * q * p * T(0.5)
        ac[i] = af - bd[i]
        bc[i] = bf + ad[i]
    end
    return (ac, bc, ad, bd)
end

@inline _coeffs_at(c::NTuple{N, AbstractVector}, n) where N = ntuple(j -> @inbounds(c[j][n]), Val(N))

"""
Indices of the four points of the stencil for output index `l` (from 1 to `N`) and
integer shift `k` (from 0 to `N - 1`), wrapped around periodically.
"""
@inline function _taps(l, k, N)
    i0 = l + k
    i0 = ifelse(i0 > N, i0 - N, i0)
    im = ifelse(i0 == 1, N, i0 - 1)
    i1 = ifelse(i0 == N, 1, i0 + 1)
    i2 = ifelse(i1 == N, 1, i1 + 1)
    return im, i0, i1, i2
end

@inline _interp(::LinearInterp, c, ym, y0, y1, y2) = y0 + c[1] * (y1 - y0)

@inline _interp(::CubicInterp, c, ym, y0, y1, y2) = c[1]*y0 + c[2]*y1 - c[3]*ym + c[4]*y2

@inline _interp(::SteffenInterp, c, ym, y0, y1, y2) = _steffen(c[1], ym, y0, y1, y2)

@inline _interp(::Weno4Interp, c, ym, y0, y1, y2) = _weno4(c[1], c[2], ym, y0, y1, y2)

"""
Monotonic cubic Hermite interpolation of Steffen (1990) on a uniform grid.
The derivatives at the two interval ends are limited so that the interpolant is monotonic
in the interval, and the final clamp guards against rounding errors, so that the
result is always between `y0` and `y1`.
"""
@inline function _steffen(t::T, ym::T, y0::T, y1::T, y2::T) where T
    sm = y0 - ym
    s0 = y1 - y0
    s1 = y2 - y1
    d0 = (sign(sm) + sign(s0)) * min(abs(sm), abs(s0), T(0.25) * abs(sm + s0))
    d1 = (sign(s0) + sign(s1)) * min(abs(s0), abs(s1), T(0.25) * abs(s0 + s1))
    result = y0 + t * (d0 + t * (3 * s0 - 2 * d0 - d1 + t * (d0 + d1 - 2 * s0)))
    return clamp(result, min(y0, y1), max(y0, y1))
end

"""
Fourth-order WENO interpolation of [Janett et al. (2019)](https://ui.adsabs.harvard.edu/abs/2019A%26A...624A.104J/abstract)
for a uniform grid, between `y0` and `y1`, a fraction `t` of the way, with `c3 = t(t - 1)/2`.
It combines the quadratic interpolants of the stencils (`ym`, `y0`, `y1`) and (`y0`, `y1`, `y2`),
with weights that favour the smoother one. The smoothness indicators are normalised by the
magnitude of the data (a scale-invariant form), avoiding overflows in single precision.
"""
@inline function _weno4(t::T, c3::T, ym::T, y0::T, y1::T, y2::T) where T
    a = y0 - ym
    b = y1 - y0
    c = y2 - y1
    q2 = y0 + t * (T(0.5) * (a + b) + T(0.5) * t * (b - a))
    q3 = y0 + t * b + c3 * (c - b)
    r = inv(abs(ym) + abs(y0) + abs(y1) + abs(y2) + floatmin(T))
    D0 = abs(c - a)
    u2 = (D0 - abs(4b - 3a - c)) * r
    u3 = (abs(a - 4b + 3c) - D0) * r
    β2 = u2 * u2 + T(WENO4_ε)
    β3 = u3 * u3 + T(WENO4_ε)
    w2 = (2 - t) * β3   # linear weights (2 - t)/3 and (1 + t)/3, divided by β
    w3 = (1 + t) * β2
    return q2 + (w3 / (w2 + w3)) * (q3 - q2)
end

"""
Shifts data along the last axis of `data_in`, writing into `data_out`, with a different
(integer, coefficients) shift for each index in the first axis, constant in the second axis.
Heights are processed in blocks, going through the shifted axis, so that all points in
a block use the same cache lines, which are then fully used. A sliding window of
four points along the shifted axis is kept for each height, so that only one new point is
loaded for each output. The windows live in work arrays with a (padded) slice for each
chunk of the second axis (one per thread), allocated once so that the loops do not allocate.
"""
function _shift_axis!(data_out, data_in, scheme, k, coeffs)
    nz, na, nl = size(data_in)
    T = eltype(data_in)
    B = SHIFT_BLOCK
    nchunks = min(Threads.nthreads(), na)
    chunk_size = cld(na, nchunks)
    windows = Array{T, 3}(undef, B + SHIFT_PAD, 4, nchunks)
    indices = Array{Int, 2}(undef, B + SHIFT_PAD, nchunks)
    Threads.@threads for c in 1:nchunks
        wm = view(windows, :, 1, c)
        w0 = view(windows, :, 2, c)
        w1 = view(windows, :, 3, c)
        w2 = view(windows, :, 4, c)
        ip = view(indices, :, c)
        @inbounds for a in ((c - 1) * chunk_size + 1):min(c * chunk_size, na)
            for n1 in 1:B:nz
                nn = min(B, nz - n1 + 1)
                for j in 1:nn
                    im, i0, i1, i2 = _taps(1, k[n1 + j - 1], nl)
                    wm[j] = data_in[n1 + j - 1, a, im]
                    w0[j] = data_in[n1 + j - 1, a, i0]
                    w1[j] = data_in[n1 + j - 1, a, i1]
                    w2[j] = data_in[n1 + j - 1, a, i2]
                    ip[j] = i2
                end
                for l in 1:nl
                    @simd for j in 1:nn
                        n = n1 + j - 1
                        data_out[n, a, l] = _interp(scheme, _coeffs_at(coeffs, n),
                                                    wm[j], w0[j], w1[j], w2[j])
                    end
                    l == nl && break
                    @simd for j in 1:nn
                        wm[j] = w0[j]
                        w0[j] = w1[j]
                        w1[j] = w2[j]
                        ip[j] = ifelse(ip[j] == nl, 1, ip[j] + 1)
                    end
                    for j in 1:nn
                        w2[j] = data_in[n1 + j - 1, a, ip[j]]
                    end
                end
            end
        end
    end
    return nothing
end
