using Interpolations

# Deterministic pseudo-random numbers in [0, 1)
pseudo(n, seed) = [mod(sin(12.9898 * seed + 78.233 * i) * 43758.5453, 1.0) for i in 1:n]

@testset "incline.jl" begin
    @testset "incline_data!" begin
        # Simple case with 4x3x3 atmosphere, z=0 at last index,
        tmp_in = zeros(Float32, 4, 3, 3)
        tmp_in[:, 2, 2] .= 1f0
        tmp_out = zeros(Float32, 4, 3, 3)
        z = collect(reverse(0:3f0))

        @test_throws ArgumentError incline_data!(tmp_in, tmp_out, z, 1., 1., 0.4, 0;
                                                 interpolation=:invalid)
        @test_throws AssertionError incline_data!(tmp_in, tmp_out, 3:-1:1, 1., 1., 0.4, 0)
        @test_throws AssertionError incline_data!(tmp_in, tmp_out, z, 0., 1., 0.4, 0)
        @test_throws AssertionError incline_data!(tmp_in, tmp_out, z, 1., 0., 0.4, 0)
        @test_throws AssertionError incline_data!(tmp_in, tmp_out, z, 1., 1., 0., 0)
        # no inclination:
        incline_data!(tmp_in, tmp_out, z, 1., 1., 1., 0)
        @test tmp_out == tmp_in

        for interp in [:linear, :cubic, :steffen, :weno4]
            # θ=π/4, ϕ=0  (along x axis), one pixel shift at each height
            incline_data!(tmp_in, tmp_out, z, 1., 1., cos(π/4), 0; interpolation=interp)
            @test tmp_out[4, :, :] == tmp_in[4, :, :]  # z=0 layer should not shift
            @test tmp_out[3, :, :] == [0 0 0
                                       0 0 1
                                       0 0 0]
            @test tmp_out[2, :, :] == [0 0 0
                                       1 0 0
                                       0 0 0]
            @test tmp_out[1, :, :] == [0 0 0
                                       0 1 0
                                       0 0 0]
            # θ=π/4, ϕ=π/2 (along y axis), one pixel shift at each height
            incline_data!(tmp_in, tmp_out, z, 1., 1., cos(π/4), π/2; interpolation=interp)
            @test tmp_out[4, :, :] == tmp_in[4, :, :]  # z=0 layer should not shift
            @test tmp_out[3, :, :] == [0 0 0
                                       0 0 0
                                       0 1 0]
            @test tmp_out[2, :, :] == [0 1 0
                                       0 0 0
                                       0 0 0]
            @test tmp_out[1, :, :] == [0 0 0
                                       0 1 0
                                       0 0 0]
            # ϕ=π/4, θ such that one pixel shift diagonally at each height
            μ = cos(atan(sqrt(2)))
            incline_data!(tmp_in, tmp_out, z, 1., 1., μ, π/4; interpolation=interp)
            @test tmp_out[4, :, :] == tmp_in[4, :, :]  # z=0 layer should not shift
            @test tmp_out[3, :, :] ≈ [0 0 0
                                      0 0 0
                                      0 0 1]
            @test tmp_out[2, :, :] ≈ [1 0 0
                                      0 0 0
                                      0 0 0]
            @test tmp_out[1, :, :] ≈ [0 0 0
                                      0 1 0
                                      0 0 0]        end
    end

    @testset "interpolation schemes" begin
        # Interpolate between y0 and y1, a fraction t of the way
        function interp(name, t, ym, y0, y1, y2)
            scheme = Muspel._interp_scheme(name)
            c = Muspel._coeffs(scheme, [t])
            return Muspel._interp(scheme, ntuple(j -> c[j][1], length(c)), ym, y0, y1, y2)
        end
        schemes = [:linear, :cubic, :steffen, :weno4]

        @test_throws "not recognised" Muspel._interp_scheme(:invalid)

        @testset "nodes and constants" begin
            for name in schemes, T in [Float32, Float64]
                y = T.(pseudo(4, 1) .+ 1)
                @test interp(name, zero(T), y...) == y[2]
                @test interp(name, one(T), y...) ≈ y[3]
                for t in T[0, 0.2, 0.5, 0.99]
                    if name == :cubic   # coefficients add up to one except for rounding
                        @test interp(name, t, y[1], y[1], y[1], y[1]) ≈ y[1]
                    else
                        @test interp(name, t, y[1], y[1], y[1], y[1]) == y[1]
                    end
                end
            end
        end

        @testset "weno4 reference values" begin
            # Values from WENO4.jl for the grid (-1, 0, 1, 2) and the stencils below, which
            # have unit sum of absolute values so that the absolute ε of WENO4.jl is the same
            # as the (relative) ε used here. Columns are t = 0, 0.25, 0.5, 0.9
            stencils = [[0.1, 0.2, 0.3, 0.4],
                        [0.0, 0.0, 0.5, 0.5],
                        [0.05, 0.4, 0.5, 0.05],
                        [0.3, -0.1, 0.2, -0.4],
                        [1.0, 0.0, 0.0, 0.0]]
            ref = [[0.2, 0.225, 0.25, 0.29],
                   [0.0, 0.1171875, 0.25, 0.456],
                   [0.4, 0.46015624999999993, 0.5, 0.5098],
                   [-0.1, -0.03759469594483811, 0.049225665185214026, 0.179498165665348],
                   [0.0, -3.281248022966831e-8, -3.1249984296109225e-8, -6.513155307310348e-9]]
            for (s, r) in zip(stencils, ref)
                @test [interp(:weno4, t, s...) for t in [0.0, 0.25, 0.5, 0.9]] ≈ r atol=1e-14
            end
            # reproduces quadratics, and is invariant to the units of the data
            p(x) = 2.5 - 1.5*x + 0.7*x^2
            for t in [0.1, 0.5, 0.8]
                @test interp(:weno4, t, p(-1), p(0), p(1), p(2)) ≈ p(t)
                @test interp(:weno4, t, 1e30 .* stencils[2]...) ≈ 1e30 * interp(:weno4, t, stencils[2]...)
            end
            # in single precision large values do not overflow
            @test isfinite(interp(:weno4, 0.5f0, 1f23, 2f23, 5f23, 6f23))
        end

        @testset "steffen" begin
            # Same as the monotonic interpolation in Interpolations.jl
            n = 30
            y = cumsum(pseudo(n, 2) .- 0.5) .* 3
            x = collect(-4.0:n+3.0)
            yp = [y[mod1(i, n)] for i in -3:n+4]
            itp = interpolate(x, yp, SteffenMonotonicInterpolation())
            for l in 2:n-2, t in [0.1, 0.5, 0.75]
                @test interp(:steffen, t, y[l-1], y[l], y[l+1], y[l+2]) ≈ itp(l - 1 + t)
            end
            # Never leaves the range of the two nearest points, whatever the contrast
            for i in 1:2000
                y = Float32.(10 .^ (20 .* pseudo(4, i) .- 10))
                t = Float32(pseudo(1, -i)[1])
                r = interp(:steffen, t, y...)
                @test min(y[2], y[3]) <= r <= max(y[2], y[3])
                @test r > 0
            end
        end

        @testset "overshoot" begin
            ts = 0.0:0.01:1.0
            # intervals next to a jump
            for step in ([0.0, 0.0, 0.0, 1.0], [0.0, 1.0, 1.0, 1.0])
                @test all(0 <= interp(:linear, t, step...) <= 1 for t in ts)
                @test minimum(interp(:cubic, t, step...) for t in ts) < -0.05 ||
                      maximum(interp(:cubic, t, step...) for t in ts) > 1.05
                @test all(0 <= interp(:steffen, t, step...) <= 1 for t in ts)
                @test all(-1e-5 < interp(:weno4, t, step...) < 1 + 1e-5 for t in ts)
            end
            # sharp drop to a small positive value
            drop = [1e8, 1.0, 1.0, 1.0]
            @test minimum(interp(:cubic, t, drop...) for t in ts) < 0
            @test minimum(interp(:steffen, t, drop...) for t in ts) == 1
            @test minimum(interp(:weno4, t, drop...) for t in ts) > -1e-6 * 1e8
        end

        @testset "accuracy" begin
            # root mean square error for a smooth function, resolved with at least 10 points
            # per wavelength (Steffen is only second order at the extrema)
            f(x) = sin(2π * x / 16 + 0.3) + 0.3*cos(2π * x / 10)
            function rms(name)
                e = Float64[]
                for l in 0:15, t in 0.0:0.1:1.0
                    push!(e, interp(name, t, f.(l-1:l+2)...) - f(l + t))
                end
                return sqrt(sum(abs2, e) / length(e))
            end
            @test rms(:cubic) < rms(:linear) / 10
            @test rms(:weno4) < rms(:linear) / 10
            @test rms(:steffen) < rms(:linear) / 2
        end
    end

    @testset "incline_data! schemes" begin
        # periodic data with a jump, in a cube that is not a multiple of the block size
        nz, ny, nx = 70, 9, 11
        z = collect(range(1000f0, -300f0, length=nz))
        data = Float32.(reshape(pseudo(nz*ny*nx, 3), nz, ny, nx) .+ 1)
        out = similar(data)
        # Direct implementation with explicit periodic indexing, to check the blocked
        # and sliding window shifting, for each axis in turn (x, then y)
        function reference(data, scheme, plan)
            nz, ny, nx = size(data)
            tmp = copy(data)
            if plan.shift_x
                for l in 1:nx, m in 1:ny, n in 1:nz
                    i0 = mod1(l + plan.kx[n], nx)
                    tmp[n, m, l] = Muspel._interp(
                        scheme, Muspel._coeffs_at(plan.cx, n), data[n, m, mod1(i0 - 1, nx)],
                        data[n, m, i0], data[n, m, mod1(i0 + 1, nx)], data[n, m, mod1(i0 + 2, nx)]
                    )
                end
            end
            result = copy(tmp)
            if plan.shift_y
                for l in 1:nx, m in 1:ny, n in 1:nz
                    i0 = mod1(m + plan.ky[n], ny)
                    result[n, m, l] = Muspel._interp(
                        scheme, Muspel._coeffs_at(plan.cy, n), tmp[n, mod1(i0 - 1, ny), l],
                        tmp[n, i0, l], tmp[n, mod1(i0 + 1, ny), l], tmp[n, mod1(i0 + 2, ny), l]
                    )
                end
            end
            return result
        end
        for interp in [:linear, :cubic, :steffen, :weno4], (μ, ϕ) in [(0.3, 0.0), (0.6, π/2), (0.4, 1.0), (0.8, -2.0)]
            incline_data!(data, out, z, 100f0, 100f0, μ, ϕ; interpolation=interp)
            scheme = Muspel._interp_scheme(interp)
            plan = Muspel._incline_plan(scheme, Float32, nz, ny, nx, z, 100f0, 100f0, μ, ϕ)
            @test out == reference(data, scheme, plan)
            buf = similar(data)
            out2 = similar(data)
            incline_data!(data, out2, z, 100f0, 100f0, μ, ϕ; interpolation=interp, buffer=buf)
            @test out2 == out
        end
        # a shift by half a pixel is the average for linear interpolation
        z2 = collect(0f0:2f0:6f0)  # shift of 0.5 pixels per height step
        d = Float32.(reshape(pseudo(60, 4), 4, 3, 5))
        o = similar(d)
        incline_data!(d, o, z2, 4f0, 1f0, cos(atan(1f0)), 0; interpolation=:linear)
        @test o[1, :, :] ≈ d[1, :, :]
        @test o[2, :, :] ≈ (d[2, :, :] .+ circshift(d[2, :, :], (0, 1))) ./ 2 rtol=1e-6
        # no overshoot when inclining a sharp front
        front = zeros(Float32, nz, ny, nx)
        front[:, :, 4:8] .= 1
        for interp in [:linear, :steffen]
            incline_data!(front, out, z, 100f0, 100f0, 0.3, 0.7; interpolation=interp)
            @test extrema(out) == (0, 1)
        end
        incline_data!(front, out, z, 100f0, 100f0, 0.3, 0.7; interpolation=:cubic)
        @test maximum(out) > 1.01
        incline_data!(front, out, z, 100f0, 100f0, 0.3, 0.7; interpolation=:weno4)
        @test -1e-4 < minimum(out) && maximum(out) < 1 + 1e-4
    end

    @testset "incline_data! allocations" begin
        # The number of allocations must not grow with the size of the data (no work arrays
        # allocated inside the loops), apart from small differences from threading
        function allocations(n, interp)
            nz = 40
            z = collect(range(1000f0, -300f0, length=nz))
            d = Float32.(reshape(pseudo(nz*n*n, 1), nz, n, n))
            o = similar(d)
            b = similar(d)
            incline_data!(d, o, z, 100f0, 100f0, 0.4, 1.0; interpolation=interp, buffer=b)
            return @allocations incline_data!(d, o, z, 100f0, 100f0, 0.4, 1.0;
                                              interpolation=interp, buffer=b)
        end
        for interp in [:linear, :cubic, :steffen, :weno4]
            @test allocations(96, interp) <= allocations(48, interp) + 10
        end
    end

    @testset "incline_atmos" begin
        FALC_multi3d = joinpath(@__DIR__, "..", "data", "atmospheres", "atm3d.FALC.3x3x82")
        FALC_multi3d_mesh = joinpath(@__DIR__, "..", "data", "atmospheres", "mesh.FALC.3x3x82")
        atm = read_atmos_multi3d(FALC_multi3d_mesh, FALC_multi3d)
        # rotation of 1D atmosphere should be the same, except z scale
        for interp in [:linear, :cubic, :steffen, :weno4]
            inc_atm = incline_atmos(atm, 0.5, 0; interpolation=interp, vector_interpolation=interp)
            same = interp == :cubic ? (a, b) -> a ≈ b : (a, b) -> a == b  # cubic has rounding
            @test same(inc_atm.electron_density, atm.electron_density)
            @test same(inc_atm.proton_density, atm.proton_density)
            @test same(inc_atm.temperature, atm.temperature)
            @test same(inc_atm.velocity.z, atm.velocity.z)
            @test all(inc_atm.z .== atm.z / 0.5)
        end
        @test_throws "not recognised" incline_atmos(atm, 0.5, 0; interpolation=:invalid)
        @test_throws "not recognised" incline_atmos(atm, 0.5, 0; vector_interpolation=:invalid)
    end

    @testset "incline_atmos!" begin
        nz, ny, nx = 6, 8, 10
        z = collect(range(2., -3., length=nz))
        x = collect(0.:nx-1) .* 2
        y = collect(0.:ny-1) .* 2
        field(seed) = reshape(pseudo(nz*ny*nx, seed), nz, ny, nx)
        scalar(seed) = 10 .^ (6 .* field(seed))   # large dynamic range
        v(seed) = (x = field(seed) .- 0.5, y = field(seed + 1) .- 0.5, z = field(seed + 2) .- 0.5)
        atm = Atmosphere(nx, ny, nz, x, y, z, scalar(5), v(6), v(9), scalar(12), scalar(13),
                         scalar(14))
        dest = deepcopy(atm)
        incline_atmos!(dest, atm, 0.6, 1.0)
        inc = incline_atmos(atm, 0.6, 1.0)
        for name in (:temperature, :electron_density, :hydrogen1_density, :proton_density)
            @test getfield(dest, name) == getfield(inc, name)
            # default scalar interpolation cannot create new extrema
            @test minimum(getfield(inc, name)) >= minimum(getfield(atm, name))
            @test maximum(getfield(inc, name)) <= maximum(getfield(atm, name))
        end
        @test dest.velocity == inc.velocity
        @test dest.magnetic_field == inc.magnetic_field
        @test dest.z == inc.z && dest.x == inc.x && dest.y == inc.y
        # input is not modified
        @test atm.z == z && atm.x == x && atm.y == y
        # cubic interpolation of scalars does create negative values
        cubic = incline_atmos(atm, 0.6, 1.0; interpolation=:cubic)
        @test minimum(cubic.temperature) < minimum(atm.temperature)
        # arrays cannot be shared
        @test_throws "must not share arrays" incline_atmos!(atm, atm, 0.6, 1.0)
        @test_throws "not recognised" incline_atmos!(deepcopy(atm), atm, 0.6, 1.0; interpolation=:invalid)
    end

    @testset "incline_atmos with magnetic field" begin
        nz, ny, nx = 4, 3, 3
        z = collect(0.:nz-1)
        x = collect(0.:nx-1)
        y = collect(0.:ny-1)
        scalar = ones(nz, ny, nx)
        v0 = (x = zeros(nz, ny, nx), y = zeros(nz, ny, nx), z = zeros(nz, ny, nx))
        b0 = (x = ones(nz, ny, nx), y = zeros(nz, ny, nx), z = zeros(nz, ny, nx))
        atm = Atmosphere(nx, ny, nz, x, y, z, scalar, v0, b0,
                         copy(scalar), copy(scalar), copy(scalar))
        @test has_magnetic_field(atm)
        # μ=1: no polar shift, φ=π/2 rotates vectors azimuthally
        inc_atm = incline_atmos(atm, 1.0, π/2)
        @test has_magnetic_field(inc_atm)
        @test all(isapprox.(inc_atm.magnetic_field.x, 0, atol=1e-15))
        @test all(inc_atm.magnetic_field.y .≈ -1)
        @test all(isapprox.(inc_atm.magnetic_field.z, 0, atol=1e-15))
        # atmospheres without 3 velocity components cannot be inclined
        atm_z = Atmosphere(nx, ny, nz, x, y, z, scalar,
                           (z = zeros(nz, ny, nx),), nothing,
                           copy(scalar), copy(scalar), copy(scalar))
        @test_throws ArgumentError incline_atmos(atm_z, 0.5, 0)
    end

    @testset "project_vector" begin
        # Azimuthal only projections
        vx = [1., 0, 0]
        vy = [0, 1., 0]
        vz = [0, 0, 1.]
        project_vector!(vx, vy, vz, 1., π/2)
        @test [vx vy vz]' ≈ [ 0  1  0
                             -1  0  0
                              0  0  1]
        vx = [1., 0, 0]
        vy = [0, 1., 0]
        vz = [0, 0, 1.]
        project_vector!(vx, vy, vz, 1., -π/2)
        @test [vx vy vz]' ≈ [ 0 -1  0
                              1  0  0
                              0  0  1]
        vx = [sqrt(2), 0]
        vy = [0, sqrt(2)]
        vz = [0, 0.]
        project_vector!(vx, vy, vz, 1., π/4)
        @test [vx vy vz]' ≈ [ 1  1
                             -1  1
                              0  0]
        vx = [sqrt(2), 0]
        vy = [0, sqrt(2)]
        vz = [0, 0.]
        project_vector!(vx, vy, vz, 1., -π/4)
        @test [vx vy vz]' ≈ [ 1 -1
                              1  1
                              0  0]
        # Polar only
        vx = [1., 0, 0]
        vy = [0, 1., 0]
        vz = [0, 0, 1.]
        project_vector!(vx, vy, vz, 0., 0)  # θ=π/2
        @test [vx vy vz]' ≈ [ 0  0 -1
                              0  1  0
                              1  0  0]
        # Polar and azimuthal
        vx = [1., 0, 0]
        vy = [0, 1., 0]
        vz = [0, 0, 1.]
        project_vector!(vx, vy, vz, 0., π)
        @test [vx vy vz]' ≈ [ 0  0 -1
                              0 -1  0
                             -1  0  0]
        vx = [0.]
        vy = [0.]
        vz = [1.]
        project_vector!(vx, vy, vz, 0.5, 0)
        @test [vx vy vz] ≈ [-sin(π/3) 0 0.5]
    end
end
