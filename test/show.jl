@testset "show.jl" begin
    @testset "Atmospheres" begin
        FALC_RH_file = joinpath(@__DIR__, "..", "data", "atmospheres", "FALC.hdf5")
        FALC_RH_nHtot_file = joinpath(@__DIR__, "..", "data", "atmospheres", "FALC_nHtot.hdf5")
        FALC_multi3d = joinpath(@__DIR__, "..", "data", "atmospheres", "atm3d.FALC.3x3x82")
        FALC_multi3d_mesh = joinpath(@__DIR__, "..", "data", "atmospheres", "mesh.FALC.3x3x82")

        atm1d = read_atmos_rh(FALC_RH_file)
        txt = sprint(show, MIME("text/plain"), atm1d)
        @test occursin("Atmosphere1D{Float32}", txt)
        @test occursin("82×3×3", txt)
        @test occursin("dims", txt)
        @test occursin("loaded in memory", txt)
        @test occursin("data size", txt)
        txt1 = sprint(show, MIME("text/plain"), atm1d[:, 1, 1])
        @test occursin("Atmosphere1D{Float32}", txt1)
        @test occursin("point", txt1)

        atm3d = read_atmos_multi3d(FALC_multi3d_mesh, FALC_multi3d)
        txt = sprint(show, MIME("text/plain"), atm3d)
        @test occursin("Atmosphere3D{Float32}", txt)
        @test occursin("3 points:", txt)
        @test occursin("dims", txt)
        @test occursin("loaded in memory", txt)
        @test occursin("data size", txt)
        txt1 = sprint(show, MIME("text/plain"), atm3d[:, 1, 1])
        @test occursin("Atmosphere1D{Float32}", txt1)
        @test occursin("point", txt1)
    end

    @testset "Atoms" begin
        H = read_atom("test_atoms/H_test.yaml")
        txt = sprint(show, MIME("text/plain"), H)
        @test occursin("AtomicModel{Float64, Int64}", txt)
        @test occursin("H (I, II)", txt)
        @test occursin("dims", txt)
        @test occursin("6×10×5", txt)
        @test occursin("levels", txt)
        @test occursin("lines", txt)
        @test occursin("continua", txt)

        txt = sprint(show, MIME("text/plain"), H.lines)
        @test occursin("Vector{AtomicLine}", txt)
        @test occursin("432", txt)
        @test occursin("λ₀", txt)

        txt = sprint(show, MIME("text/plain"), H.lines[5])
        @test occursin("AtomicLine{Float64}", txt)
        @test occursin("656.5", txt)
        @test occursin("1 Zeeman component", txt)
        @test occursin("metadata", txt)

        txt = sprint(show, MIME("text/plain"), H.continua)
        @test occursin("Vector{AtomicContinuum}", txt)
        @test occursin("100", txt)
        @test occursin("⩘ λ", txt)

        txt = sprint(show, MIME("text/plain"), H.continua[1])
        @test occursin("AtomicContinuum{Float64}", txt)
        @test occursin("91.18", txt)
        @test occursin("metadata", txt)
    end
end
