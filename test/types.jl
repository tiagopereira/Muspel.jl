@testset "types.jl" begin
    @testset "Atmosphere" begin
        x64 = ones(Float64, 2)
        z64 = collect(0.:1.:6.)
        temp64 = ones(Float64, 7, 2, 2) # (z,y,x)
        h64 = ones(Float64, 7, 2, 2, 3) # (z,y,x,h)
        x32 = convert.(Float32, x64)
        z32 = convert.(Float32, z64)
        temp32 = convert.(Float32, temp64)
        h32 = convert.(Float32, h64)
        # Input types:
        @test_throws MethodError Atmosphere(x32, x64, z64, temp64, temp64, temp64, h64)
        # Input dimensions:
        @test_throws AssertionError Atmosphere(x64, x64, x64, temp64, temp64, temp64, h64)
        @test_throws AssertionError Atmosphere(
                                            x64, x64, z64, ones(1,1,1), temp64, temp64, h64)
        @test_throws AssertionError Atmosphere(
                                            x64, x64, z64, temp64, ones(1,1,1), temp64, h64)
        @test_throws AssertionError Atmosphere(
                                            x64, x64, z64, temp64, temp64, ones(1,1,1), h64)
        # Single precision/Double precision
        y64 = x64 .+ 1
        vz64 = (temp64 .+ 2)
        ne64 = (temp64 .+ 3)
        atm64 = Atmosphere(x64, y64, z64, temp64, vz64, ne64, h64)
        atm32 = Atmosphere(x32, x32, z32, temp32, temp32, temp32, h32)
        fields = fieldnames(Atmosphere)
        sysInt = typeof(1)
        ftypes64 = [eltype(getfield(atm64, f)) for f in fields]
        @test all( ((ftypes64 .== sysInt) + (ftypes64 .== Float64)) .== 1)
        ftypes32 = [eltype(getfield(atm32, f)) for f in fieldnames(Atmosphere)]
        @test all( ((ftypes32 .== sysInt) + (ftypes32 .== Float32)) .== 1)
        # Output fields:
        @test atm64.nx == length(x64)
        @test atm64.ny == length(y64)
        @test atm64.nz == length(z64)
        @test atm64.nh_levels == size(h64)[end]
        @test atm64.x == x64
        @test atm64.y == y64
        @test atm64.z == z64
        @test atm64.temperature == temp64
        @test atm64.velocity_z == vz64
        @test atm64.electron_density == ne64
        @test atm64.hydrogen_density == h64
    end
    @testset "AtomicLine" begin
        @test fieldnames(AtomicLine) == (:nλ, :χup, :χlo, :gup, :glo, :Aul, :Blu, :Bul, :λ0,
                                         :f_value, :λ, :PRD, :Voigt, :label_up, :label_lo,
                                         :γ)
    end
    @testset "AtomicContinuum" begin
        @test fieldnames(AtomicContinuum) == (:up, :lo, :nλ, :λedge, :σ, :λ)
    end
    @testset "AtomicModel" begin
        # Fieldnames unchanged (will have to add collisions)
        @test fieldnames(AtomicModel) == (:element, :nlevels, :nlines, :ncontinua, :Z,
                                          :mass, :χ, :g, :stage, :label, :lines, :continua)
    end
end
