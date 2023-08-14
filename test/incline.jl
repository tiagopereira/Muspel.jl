@testset "incline.jl" begin
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
