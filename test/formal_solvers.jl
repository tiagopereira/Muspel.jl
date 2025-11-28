@testset "formal_solvers.jl" begin
    @testset "Weights" begin
        @test all(Muspel._w3(60.) .≈ (0, 0.016666666666666666, 0.98333333333333))
        @test all(Muspel._w3(1.) .≈ (exp(-1), 1 - 2*exp(-1), exp(-1)))
        @test all(Muspel._w3(1f-6) .≈ (0.999999f0, 4.9999966f-7, 4.999998f-7))
        @test sum(Muspel._w3(50.)) == 1
        @test sum(Muspel._w3(0.1)) == 1
        @test sum(Muspel._w3(1f-4)) == 1
    end
    @testset "Piecewise" begin
        # Constant source function
        z = collect(LinRange(1, 1e6, 20))
        alpha = ones(20) * 1e-20
        S = ones(20)
        @test piecewise_1D_linear(z, alpha, S) ≈ S
        @test begin
            result = similar(S)
            piecewise_1D_linear!(z, alpha, S, result)
            result ≈ S
        end
        @test begin
            result = similar(S)
            piecewise_1D_bezier3!(z, alpha, S, result)
            result ≈ S
        end
        alpha = ones(20)
        @test piecewise_1D_linear(z, alpha, S) ≈ S
        @test piecewise_1D_linear(z, alpha, S;
                                  initial_condition=:zero)[[1, end]] ≈ [S[1], S[1]*0]
        @test begin
            result = similar(S)
            piecewise_1D_linear!(z, alpha, S, result)
            result ≈ S
        end
        @test begin
            result = similar(S)
            piecewise_1D_bezier3!(z, alpha, S, result)
            result ≈ S
        end
        @test begin
            result = similar(S)
            piecewise_1D_linear!(z, alpha, S, result; initial_condition=:zero)
            result[[1, end]] ≈ [S[1], S[1]*0]
        end
        @test begin
            result = similar(S)
            piecewise_1D_bezier3!(z, alpha, S, result; initial_condition=:zero)
            result[[1, end]] ≈ [S[1], S[1]*0]
        end

        # Linear extinction and source function, test reversibility
        alpha = collect(LinRange(1e-3, 1e-5, 20))
        S = collect(LinRange(1, 100, 20))
        @test (piecewise_1D_linear(z, reverse(alpha), reverse(S); to_end=true) ≈
               reverse(piecewise_1D_linear(z, alpha, S)))
        @test begin
            result_i = similar(S)
            piecewise_1D_linear!(z, reverse(alpha), reverse(S), result_i; to_end=true)
            result = similar(S)
            piecewise_1D_linear!(z, alpha, S, result)
            result_i ≈ reverse(result)
        end
        @test begin
            result_i = similar(S)
            piecewise_1D_bezier3!(z, reverse(alpha), reverse(S), result_i; to_end=true)
            result = similar(S)
            piecewise_1D_bezier3!(z, alpha, S, result)
            result_i ≈ reverse(result)
        end

        # Exceptions
        result = similar(S)
        @test_throws ErrorException piecewise_1D_linear(z, alpha, S; initial_condition=:aaa)
        @test_throws ErrorException piecewise_1D_linear!(z, alpha, S, result; initial_condition=:aaa)
        @test_throws ErrorException piecewise_1D_bezier3!(z, alpha, S, result; initial_condition=:aaa)
    end
    @testset "Feautrier" begin
        z = collect(LinRange(2e6, -1e5, 20))
        alpha = 1e-5 * ones(20)
        S = zeros(20)
        # Simple tests
        @test feautrier(z, alpha, S) ≈ S
        S = 100 * ones(20)
        # Against implementation
        @test feautrier(z, alpha, S)[1] * 2 ≈ 106.65292755967045
        # Test reversibility
        @test feautrier(z, alpha, S)[1] ≈ feautrier(z, reverse(alpha), reverse(S))[end]
    end
    @testset "Bezier3" begin
        # Testing against reference C implementation
        res = [a for a in Muspel.bezier3_coeffs(0.01)]
        @test res ≈ [
            0.0024800833333333333,
            0.0024950083333333335,
            0.0024850500000000004,
            0.002490025,
            0.9900498333333333
        ]
        res = [a for a in Muspel.bezier3_coeffs(1.)]
        @test res ≈ [
            0.11392894125692266,
            0.207276647028654,
            0.1400215586575957,
            0.17089341188538398,
            0.36787944117144233
        ]
        res = [a for a in Muspel.bezier3_coeffs(100.)]
        @test res ≈ [6e-06, 0.970594, 0.000582, 0.028818, 0.0]
        # Test precision keeping
        res = [a for a in Muspel.bezier3_coeffs(1f-2)]
        @test res ≈ [
            0.0024800833333333333,
            0.0024950083333333335,
            0.0024850500000000004,
            0.002490025,
            0.9900498333333333
        ]
        res = [a for a in Muspel.bezier3_coeffs(1f0)]
        @test res ≈ [
            0.11392894125692266,
            0.207276647028654,
            0.1400215586575957,
            0.17089341188538398,
            0.36787944117144233
        ]
        res = [a for a in Muspel.bezier3_coeffs(1f2)]
        @test res ≈ [6e-06, 0.970594, 0.000582, 0.028818, 0.0]
        # Test type stability
        @test Muspel.bezier3_coeffs(1f0) isa NTuple{5, Float32}
        @test Muspel.bezier3_coeffs(1e0) isa NTuple{5, Float64}

        # Testing simple cases
        @test Muspel.cent_deriv(1, 1, 1, 1, 1) == 0.
        @test Muspel.cent_deriv(1, 1, 1, 2, 3) == 1.
        @test Muspel.cent_deriv(1, 1, 30, 20, 10) == -10.
        @test isnan(Muspel.cent_deriv(0, 1, 1, 1, 1))
        @test isnan(Muspel.cent_deriv(1, 0, 1, 1, 1))
    end
end
