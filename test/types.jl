@testset "types.jl" begin
    @testset "Atmosphere" begin
        tmp = ones(Float64, 10) * 2
        atm1 = Atmosphere1D(1, 1, 10, tmp, tmp, tmp, tmp, tmp, tmp)
        @test atm1 isa Atmosphere1D{1, Float64}
        tmp2 = repeat([1. 2. 3.], outer=[10,1,3])
        atm2 = Atmosphere1D(3, 3, 10, tmp, tmp2, tmp2, tmp2, tmp2, tmp2)
        @test atm2 isa Atmosphere1D{3, Float64}
        @test atm2[1,1] isa Atmosphere1D{1, Float64}
        @test atm2[2,2].temperature == atm1.temperature
        @test_throws MethodError atm2[1]
        atm3 = Atmosphere3D(3, 3, 10, tmp, tmp, tmp, tmp2, tmp2, tmp2, tmp2, tmp2, tmp2, tmp2)
        @test atm3[2,2].temperature == atm2[2,2].temperature
        @test atm3[3,3] isa Atmosphere1D{1, Float64}
        tmp3 = tmp2[:, :, 1]
        atm4 = Atmosphere1D(1, 3, 10, tmp, tmp3, tmp3, tmp3, tmp3, tmp3)
        @test atm4 isa Atmosphere1D{2, Float64}
        @test atm4[1] isa Atmosphere1D{1, Float64}
        @test atm4[2].proton_density == atm1.proton_density
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
