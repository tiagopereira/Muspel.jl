@testset "utils.jl" begin
    @test blackbody_λ(500, 0.0) == 0.0
    @test isnan(blackbody_λ(0, 5e3))
    @test blackbody_λ(5, 1f5) ≈ 12.134454f0
    @test blackbody_λ(100, 1f4) ≈ 6.7204614
    @test blackbody_λ(500, 6f3) ≈ 31.756907
    @test blackbody_λ(5000, 1f7) ≈ 132.43155
    @test blackbody_λ(50000, 1f2) ≈ 2.2726104f-8
end
