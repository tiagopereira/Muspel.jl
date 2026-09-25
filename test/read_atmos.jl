@testset "read_atmos.jl" begin
    @testset "RH optional datasets and timesteps" begin
        using HDF5
        mktempdir() do dir
            testfile = joinpath(dir, "rh_synth.hdf5")
            nz, ny, nx, nt = 4, 2, 3, 2
            x = Float32[1, 2, 3]
            y = Float32[4, 5]
            z = Float32[i + 10 * (t - 1) for i in 1:nz, t in 1:nt]
            temp = zeros(Float32, nz, ny, nx, nt)
            temp[:, :, :, 1] .= 5000
            temp[:, :, :, 2] .= 6000
            ne = fill(1f20, nz, ny, nx, nt)
            vx = fill(1f3, nz, ny, nx, nt); vx[:, :, :, 2] .= 2f3
            vy = fill(3f3, nz, ny, nx, nt)
            vz = fill(5f3, nz, ny, nx, nt); vz[:, :, :, 2] .= 6f3
            bx = fill(0.1f0, nz, ny, nx, nt)
            by = fill(0.2f0, nz, ny, nx, nt)
            bz = fill(0.3f0, nz, ny, nx, nt); bz[:, :, :, 2] .= 0.4f0
            htot = fill(1f25, nz, ny, nx)
            pops = Array{Float32}(undef, nz, ny, nx, 1, nt)
            pops .= 1f25
            h5write(testfile, "x", x)
            h5write(testfile, "y", y)
            h5write(testfile, "z", z)
            h5write(testfile, "temperature", temp)
            h5write(testfile, "electron_density", ne)
            h5write(testfile, "velocity_x", vx)
            h5write(testfile, "velocity_y", vy)
            h5write(testfile, "velocity_z", vz)
            h5write(testfile, "B_x", bx)
            h5write(testfile, "B_y", by)
            h5write(testfile, "B_z", bz)
            h5write(testfile, "hydrogen_populations", pops)

            # defaults: first timestep, z velocity only, no magnetic field
            atm = read_atmos_rh(testfile)
            @test keys(atm.velocity) == (:z,)
            @test !has_magnetic_field(atm)
            @test all(atm.temperature .== 5000f0)
            @test atm.z == z[:, 1]
            @test all(atm.velocity.z .== 5f3)
            @test atm.proton_density .+ atm.hydrogen1_density ≈ htot
            # select timestep, full velocity, magnetic field
            atm2 = read_atmos_rh(testfile; index=2, read_fullv=true, read_B=true)
            @test keys(atm2.velocity) == (:x, :y, :z)
            @test all(atm2.temperature .== 6000f0)
            @test atm2.z == z[:, 2]
            @test all(atm2.velocity.x .== 2f3)
            @test all(atm2.velocity.y .== 3f3)
            @test all(atm2.velocity.z .== 6f3)
            @test has_magnetic_field(atm2)
            @test all(atm2.magnetic_field.x .== 0.1f0)
            @test all(atm2.magnetic_field.y .== 0.2f0)
            @test all(atm2.magnetic_field.z .== 0.4f0)
            @test atm2.proton_density .+ atm2.hydrogen1_density ≈ htot
            # errors
            @test_throws ArgumentError read_atmos_rh(testfile; index=3)
            falc = joinpath(@__DIR__, "..", "data", "atmospheres", "FALC.hdf5")
            @test_throws ArgumentError read_atmos_rh(falc; read_B=true)
            @test_throws ArgumentError read_atmos_rh(falc; read_fullv=true)
        end
    end

    FALC_RH_file = joinpath(@__DIR__, "..", "data", "atmospheres", "FALC.hdf5")
    FALC_RH_nHtot_file = joinpath(@__DIR__, "..", "data", "atmospheres", "FALC_nHtot.hdf5")
    FALC_multi3d = joinpath(@__DIR__, "..", "data", "atmospheres", "atm3d.FALC.3x3x82")
    FALC_multi3d_mesh = joinpath(@__DIR__, "..", "data", "atmospheres", "mesh.FALC.3x3x82")
    FALC_multi3d_mesh2 = joinpath(@__DIR__, "..", "data", "atmospheres", "mesh.FALC.3x3x82_lines")

    @testset "RH" begin
        # Atmosphere with z velocity only
        atm = read_atmos_rh(FALC_RH_file)
        @test atm isa Atmosphere{3, Float32, Array{Float32, 3}, Vector{Float32}}
        @test keys(atm.velocity) == (:z,)
        @test !has_magnetic_field(atm)
        @test typeof(atm[:, 1, 1]) <: Atmosphere{1, <:AbstractFloat}
        @test typeof(atm[:, :, 1]) <: Atmosphere{2, <:AbstractFloat}
        @test typeof(atm[:, 1, :]) <: Atmosphere{2, <:AbstractFloat}
        @test typeof(atm[1:3, 1:3, 1:3]) <: Atmosphere{3, <:AbstractFloat}
        @test atm[:, 1, 1].temperature == atm.temperature[:, 1, 1]
        @test atm[:, 1, 1].velocity.z == atm.velocity.z[:, 1, 1]
        @test atm[:, 1, 1].electron_density == atm.electron_density[:, 1, 1]
        @test atm[:, 1, 1].hydrogen1_density == atm.hydrogen1_density[:, 1, 1]
        @test atm[:, 1, 1].proton_density == atm.proton_density[:, 1, 1]
        @test atm[:, 1, 1].z == atm.z
        @test maximum(atm[:, 1, 1].temperature) == 1e5
        @test minimum(atm[:, 1, 1].temperature) == 4.5f3
        # Compare explicit timestep selection
        atm2 = read_atmos_rh(FALC_RH_file; index=1)
        @test all(atm.temperature .== atm2.temperature)
        # Compare with version with nHtot
        atmH = read_atmos_rh(FALC_RH_nHtot_file)
        @test all(atmH.velocity.z .== atm.velocity.z)
        @test all(atmH.temperature .== atm.temperature)
        @test all(atmH.electron_density .== atm.electron_density)
        nH1 = atmH.proton_density[:, 1, 1] .+ atmH.hydrogen1_density[:, 1, 1]
        nH2 = atm.proton_density[:, 1, 1] .+ atm.hydrogen1_density[:, 1, 1]
        @test nH2 ≈ nH1
    end

    @testset "Multi3D" begin
        atm = read_atmos_multi3d(FALC_multi3d_mesh, FALC_multi3d)
        @test atm isa Atmosphere{3, Float32, Array{Float32, 3}, Vector{Float32}}
        @test keys(atm.velocity) == (:x, :y, :z)
        @test !has_magnetic_field(atm)
        @test atm[:, 1, 1].temperature == atm.temperature[:, 1, 1]
        @test atm[:, 1, 1].velocity.z == atm.velocity.z[:, 1, 1]
        @test atm[:, 1, 1].velocity.x == atm.velocity.x[:, 1, 1]
        @test atm[:, 1, 1].electron_density == atm.electron_density[:, 1, 1]
        @test atm[:, 1, 1].hydrogen1_density == atm.hydrogen1_density[:, 1, 1]
        @test atm[:, 1, 1].proton_density == atm.proton_density[:, 1, 1]
        @test atm[:, 1, 1].z == atm.z
        @test maximum(atm[:, 1, 1].temperature) == 1e5
        @test minimum(atm[:, 1, 1].temperature) == 4.5f3
        # Compare against FALC in RH format
        atm_RH = read_atmos_rh(FALC_RH_file)
        @test atm.velocity.z == atm_RH.velocity.z
        @test atm.temperature == atm_RH.temperature
        @test atm.electron_density ≈ atm_RH.electron_density
        atm_RH = read_atmos_rh(FALC_RH_nHtot_file)
        @test atm[:, 1, 1].proton_density ≈ atm_RH[:, 1, 1].proton_density
        @test atm[:, 1, 1].hydrogen1_density ≈ atm_RH[:, 1, 1].hydrogen1_density
        # Tests of mesh, arrays in single line vs split lines
        z1 = Muspel.read_mesh(FALC_multi3d_mesh)[end]
        z2 = Muspel.read_mesh(FALC_multi3d_mesh2)[end]
        @test z1 == z2
    end
end
