@testset "types.jl" begin
    @testset "Atmosphere" begin
        tmp = ones(Float64, 10) * 2
        atm1 = Atmosphere(1, 1, 10, Float64[], Float64[], tmp, tmp,
                          (z = tmp,), nothing, tmp, tmp, tmp)
        @test atm1 isa Atmosphere{1, Float64}
        @test keys(atm1.velocity) == (:z,)
        @test !has_magnetic_field(atm1)
        tmp2 = repeat([1. 2. 3.], outer=[10,1,3])
        atm3 = Atmosphere(3, 3, 10, Float64[], Float64[], tmp, tmp2,
                          (z = tmp2,), nothing, tmp2, tmp2, tmp2)
        @test atm3 isa Atmosphere{3, Float64}
        @test atm3[1:3, 1:3, 1:3] isa typeof(atm3)
        @test typeof(atm3[:, 1, 1]) <: Atmosphere{1, Float64}
        @test typeof(atm3[:, :, 1]) <: Atmosphere{2, Float64}
        @test typeof(atm3[:, 1, :]) <: Atmosphere{2, Float64}
        @test typeof(atm3[:, 1, :][:, 1]) <: Atmosphere{1, Float64}
        @test typeof(atm3[:, 1, :][:, 1][1:3]) <: Atmosphere{1, Float64}
        @test atm3[:, 2, 2].temperature == atm1.temperature
        @test_throws ArgumentError atm3[1]
        @test_throws ArgumentError atm3[1, :]
        @test_throws ArgumentError atm3[:, 1]
        @test_throws ArgumentError atm3[:, 1, 1, 1]
        @test_throws ArgumentError atm3[:, 1, []]

        atm3D = Atmosphere(3, 3, 10, tmp, tmp, tmp, tmp2,
                           (x = tmp2, y = tmp2, z = tmp2), nothing,
                           tmp2, tmp2, tmp2)
        @test keys(atm3D.velocity) == (:x, :y, :z)
        @test atm3D[:, 2, 2].temperature == atm3[:, 2 ,2].temperature
        @test atm3D[1:3, 1:3, 1:3] isa typeof(atm3D)
        # reduced slices lower N but keep all velocity components
        @test typeof(atm3D[:, 1, 1]) <: Atmosphere{1, Float64}
        @test typeof(atm3D[:, 1:3, 1]) <: Atmosphere{2, Float64}
        @test typeof(atm3D[:, 1, 1:3]) <: Atmosphere{2, Float64}
        @test typeof(atm3D[:, 1:3, 1]) <: Atmosphere{2, Float64}
        @test typeof(atm3D[1, 1:3, 1:3]) <: Atmosphere{2, Float64}
        @test typeof(atm3D[1, 1, 1:3]) <: Atmosphere{1, Float64}
        @test typeof(atm3D[1, 1:3, 1]) <: Atmosphere{1, Float64}
        @test keys(atm3D[:, 1, 1].velocity) == (:x, :y, :z)
        @test atm3D[:, 1, 1].velocity.x == atm3D.velocity.x[:, 1, 1]
        @test atm3D[:, 1, 1].velocity.z == atm3D.velocity.z[:, 1, 1]
        # horizontal slices rotate the axis vectors only
        h = atm3D[1, 1:3, 1:3]
        @test h.nz == 3 && h.ny == 1 && h.nx == 3
        @test h.z == atm3D.y[1:3]
        @test h.velocity.y == atm3D.velocity.y[1, 1:3, 1:3]
        @test_throws ArgumentError atm3D[:, 1, []]
        @test_throws ArgumentError atm3D[1, 1, 1]

        tmp3 = tmp2[:, :, 1]
        atm4 = Atmosphere(1, 3, 10, Float64[], tmp, tmp, tmp3,
                          (z = tmp3,), nothing, tmp3, tmp3, tmp3)
        @test atm4 isa Atmosphere{2, Float64}
        @test atm4[:, 1] isa Atmosphere{1, Float64}
        @test atm4[:, 2].proton_density == atm1.proton_density
    end
    @testset "Atmosphere with magnetic field" begin
        tmp = ones(Float64, 10) * 2
        tmp2 = repeat([1. 2. 3.], outer=[10,1,3])
        bfield = (x = tmp2, y = tmp2, z = tmp2)
        atm = Atmosphere(3, 3, 10, tmp, tmp, tmp, tmp2,
                         (x = tmp2, y = tmp2, z = tmp2), bfield,
                         tmp2, tmp2, tmp2)
        @test has_magnetic_field(atm)
        col = atm[:, 1, 1]
        @test col isa Atmosphere{1, Float64}
        @test has_magnetic_field(col)  # B survives slicing
        @test col.magnetic_field.x == tmp2[:, 1, 1]
        @test col.velocity.y == tmp2[:, 1, 1]
        # invalid component types or keys are caught by the constructor
        @test_throws ArgumentError Atmosphere(
            1, 1, 10, Float64[], Float64[], tmp, tmp,
            (z = ones(10, 2),), nothing, tmp, tmp, tmp)
        @test_throws ArgumentError Atmosphere(
            1, 1, 10, Float64[], Float64[], tmp, tmp,
            (w = tmp,), nothing, tmp, tmp, tmp)
        @test_throws ArgumentError Atmosphere(
            3, 3, 10, tmp, tmp, tmp, tmp2,
            (x = tmp2, y = tmp2, z = tmp2), (z = tmp2,),
            tmp2, tmp2, tmp2)
    end
    @testset "AtomicContinuum" begin
        @test fieldnames(AtomicContinuum) == (:up, :lo, :nλ, :λedge, :σ, :λ)
    end
    @testset "AtomicModel" begin
        # Fieldnames unchanged (will have to add collisions)
        @test fieldnames(AtomicModel) == (:element, :nlevels, :nlines, :ncontinua, :Z,
                                          :mass, :χ, :g, :stage, :label, :lines, :continua)
    end
    @testset "RTBuffer" begin
        ndep = 5
        nwave = 2
        tmp = RTBuffer(ndep, nwave, Float32)
        @test size(tmp.intensity) == (nwave,)
        @test size(tmp.source_function) == (ndep,)
    end
end
