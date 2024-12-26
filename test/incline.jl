@testset "incline.jl" begin
    @testset "incline_data!" begin
        # Simple case with 4x3x3 atmosphere, z=0 at last index,
        tmp_in = zeros(Float32, 4, 3, 3)
        tmp_in[:, 2, 2] .= 1f0
        tmp_out = zeros(Float32, 4, 3, 3)
        tmp_inv = permutedims(tmp_in, (3, 2, 1))
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

        for interp in [:linear, :cubic]
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
            tmp_inv = permutedims(tmp_in, (3, 2, 1))
            incline_data_inv!(tmp_inv, z, 1., 1., cos(π/4), 0; interpolation=interp)
            @test permutedims(tmp_inv, (3, 2, 1)) ≈ tmp_out

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
            tmp_inv = permutedims(tmp_in, (3, 2, 1))
            incline_data_inv!(tmp_inv, z, 1., 1., cos(π/4), π/2; interpolation=interp)
            @test permutedims(tmp_inv, (3, 2, 1)) == tmp_out

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
                                      0 0 0]
            tmp_inv = permutedims(tmp_in, (3, 2, 1))
            incline_data_inv!(tmp_inv, z, 1., 1., μ, π/4; interpolation=interp)
            @test permutedims(tmp_inv, (3, 2, 1)) == tmp_out
        end
    end

    @testset "incline_atmos" begin
        FALC_multi3d = joinpath(@__DIR__, "..", "data", "atmospheres", "atm3d.FALC.3x3x82")
        FALC_multi3d_mesh = joinpath(@__DIR__, "..", "data", "atmospheres", "mesh.FALC.3x3x82")
        atm = read_atmos_multi3d(FALC_multi3d_mesh, FALC_multi3d)
        # rotation of 1D atmosphere should be the same, except z scale
        inc_atm = incline_atmos(atm, 0.5, 0)
        @test all(inc_atm.electron_density .== atm.electron_density)
        @test all(inc_atm.proton_density .== atm.proton_density)
        @test all(inc_atm.temperature .== atm.temperature)
        @test all(inc_atm.velocity_z .== atm.velocity_z)
        @test all(inc_atm.z .== atm.z / 0.5)
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
