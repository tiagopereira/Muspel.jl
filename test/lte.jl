using Muspel
using Test

@testset "lte.jl" begin
    H_empty = read_atom("test_atoms/H_test_empty.yaml")
    populations = saha_boltzmann(H_empty, 1000.0, 1e15)
    @test length(populations) == H_empty.nlevels
    @test_throws MethodError saha_boltzmann(H_empty, 1, 1e15)
    @test_throws MethodError saha_boltzmann(H_empty, 1e6, 1)
    # Limits:
    @test all(saha_boltzmann(H_empty, 1e2, 1e15) .== [1.0, 0, 0, 0])
    @test all(isapprox.(saha_boltzmann(H_empty, 1e20, 1e15), [0, 0, 0, 1.0], atol=1e-20))
    @test sum(saha_boltzmann(H_empty, 1e3, 1e20)) ≈ 1
    @test sum(saha_boltzmann(H_empty, 1e6, 1e20)) ≈ 1
    # Two versions:
    @test saha_boltzmann(
        H_empty.χ, H_empty.g, H_empty.stage, 6e3, 1e18
    ) == saha_boltzmann(H_empty, 6e3, 1e18)
    # Output against previous implementation:
    temp = [6e3, 1e4, 1e6]
    ne = [1e15, 1e17, 1e20]
    @test all(saha_boltzmann(H_empty, temp[1], ne[1]) .≈ [
        0.19060553724835108,
        2.0684680211781896e-9,
        1.2062400692324476e-10,
        0.8093944605625569
    ])
    @test all(saha_boltzmann(H_empty, temp[2], ne[1]) .≈ [
        2.9537651421186625e-6,
        8.560673553408263e-11,
        2.152047770595797e-11,
        0.9999970461277307
    ])
    @test all(saha_boltzmann(H_empty, temp[3], ne[1]) .≈ [
        4.849223956870989e-16,
        1.723189157094351e-15,
        3.793124397685543e-15,
        0.9999999999999939
    ])
    @test all(saha_boltzmann(H_empty, temp[1], ne[2]) .≈ [
        0.9592653859485454,
        1.0410032170641385e-8,
        6.070675397279802e-10,
        0.040734603034355055
    ])
    @test all(saha_boltzmann(H_empty, temp[1], ne[3]) .≈ [
        0.9999575259463137,
        1.0851626845768895e-8,
        6.328194095197636e-10,
        4.246256924027464e-5
    ])
end
