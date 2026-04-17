using AtomicData
using Interpolations
using Unitful
import Muspel: create_σ_itp_LTE, create_σ_itp_NLTE, get_atoms_bf_interpolant, σH_continuum, σH_atoms_bf

@testset "background.jl" begin
    abund = get_solar_abundances()
    He = read_atom("test_atoms/He_test.yaml")
    H = read_atom("test_atoms/H_test.yaml")
    H_empty = read_atom("test_atoms/H_test_empty.yaml")
    atoms = Vector{AtomicModel}([He, H, H_empty])
    atoms_empty = Vector{AtomicModel}(undef, 0)
    σ_atom_itp = get_atoms_bf_interpolant(atoms)
    σ_atom_itp_empty = get_atoms_bf_interpolant(atoms_empty)
    λ = [92., 500., 850., 1600,  2000.]
    log_ne = 15.0:0.05:20.0
    log_temp = 3.0:0.015:6.0
    itp_lte = create_σ_itp_LTE(λ[1], log_temp, log_ne, H, atoms, σ_atom_itp)
    itp_nlte = create_σ_itp_NLTE(λ[1], log_temp, log_ne, atoms, σ_atom_itp)
    itp_nlte_empty = create_σ_itp_NLTE(λ[1], log_temp, log_ne, atoms_empty, σ_atom_itp_empty)
    temp = [2050.0, 5770.0, 10100.0]
    n = [1.0e15, 7.1e17, 1.0e20]

    @testset "σ_itp" begin
        @test isa(itp_lte, ExtinctionItpLTE{Float64})
        @test isa(itp_nlte, ExtinctionItpNLTE{Float64})
        @test itp_nlte_empty.σ_atoms.(log_temp, 19) ≈ zeros(Float64, length(log_temp))
        @test_throws MethodError create_σ_itp_LTE(λ[2], [3.1, 3.3], log_ne, H, atoms, σ_atom_itp)
        @test_throws MethodError create_σ_itp_LTE(λ[2], log_temp, [15, 19], H, atoms, σ_atom_itp)
        @test_throws MethodError create_σ_itp_NLTE(λ[2], [3.1, 3.3], log_ne, atoms, σ_atom_itp)
        @test_throws MethodError create_σ_itp_NLTE(λ[2], log_temp, [15, 19], atoms, σ_atom_itp)
        log_ne2 = 15.0:0.025:20.0 # make same number of elements as log_temp
        itp_nlte2 = create_σ_itp_NLTE(λ[2], log_temp, log_ne2, atoms, σ_atom_itp)
        atoms_list = [
            "test_atoms/He_test.yaml",
            "test_atoms/H_test.yaml",
            "test_atoms/H_test_empty.yaml"
        ]
        z_tmp = [1e5, 0., -1e5]
        atmos_test = Atmosphere1D(1, 1, length(log_temp), z_tmp, 10 .^ log_temp,
                                  [0., 0., 0.], 10 .^log_ne2, n, n)
        itp_test = get_σ_itp(atmos_test, 500., atoms_list; npts=101)
        @test isa(itp_test, ExtinctionItpNLTE{Float64})
        @test all(isapprox.(
            α_cont(itp_nlte2, 5000., 1e20, 1e20, 1e20),
            α_cont(itp_test, 5000., 1e20, 1e20, 1e20),
            rtol=1e-3
        ))
        @test_throws SystemError get_σ_itp(atmos_test, 500., ["nofile.yaml"])
    end

    @testset "α_cont" begin
        for ni in n
            ion_frac = [saha_boltzmann(H, t, ni, 1.)[end] for t in temp]
            nHI = ni .* (1 .- ion_frac)
            nHII = ni .* ion_frac
            # Check if LTE and NLTE give same result, assuming Saha for hydrogen
            @test isapprox(
                stack(α_cont.(Ref(itp_lte), temp, ni, ni)),
                stack(α_cont.(Ref(itp_nlte), temp, ni, nHI, nHII)),
                rtol=1e-5,
            )

            # Check if interpolant with no atoms gives same result as α_cont_no_itp
            @test isapprox(
                stack(α_cont_no_itp.(λ[1], temp, ni, nHI, nHII)),
                stack(α_cont.(Ref(itp_nlte_empty), temp, ni, nHI, nHII)),
                rtol=1e-5,
            )
        end
        # Some checks against implementation
        α_prev = [3.099735426728638e-9  5.737579705215677e-11  2.1702740948438375e-10
                  2.329923339815215e-8  3.531681790200384e-10  1.4127711713968792e-9
                  3.188362344689906e-8  4.5674908136823566e-10 1.9728523541442412e-9
                  1.1169536097258309e-9 1.5271949210413109e-10 5.402565280664961e-9
                  1.426592398017141e-10 2.2211042523024984e-10 9.191544796655046e-9 ]
        j_prev = 6.652458732173518e-11
        for (i, λi) in enumerate(λ)
            itp = create_σ_itp_LTE(λi, log_temp, log_ne, H, atoms, σ_atom_itp)
            res = stack(α_cont.(Ref(itp), temp, 1e18, 1e20))
            @test res[1, :] ≈ α_prev[i, :, :]
            @test all(res[2, :] .≈ j_prev)
        end
    end

    @testset "α_cont_no_itp" begin
        # Against previous implementation:
        prev = [0.3875384794861622     0.6652458732173518
                2.906477474892972      0.6662737892007118
                3.9632382586598986     0.6653610542464922
                1.307251908032521      0.66525482413674265
                1.876491725439709      0.6652495268358745] * 1e-8
        α_prev = [0.3875384794861622, 2.906477474892972, 3.9632382586598986,
                  1.307251908032521, 1.876491725439709] * 1e-8
        j_prev = [0.6652458732173518, 0.6662737892007118, 0.6653610542464922,
                  0.66525482413674265, 0.6652495268358745] * 1e-8
        res = stack(α_cont_no_itp.(λ, 6e3, 1e20, 1e20, 4e15))
        @test res[1, :] ≈ α_prev
        @test res[2, :] ≈ j_prev
        # Consistency check, Hmin maximum extinction (thermal part)
        @test argmax(res[1, :]) == 3
    end

    @testset "σH" begin
        # Against previous implementation:
        prev = [0.38752066260764906, 2.9073809748885355, 3.963188644121838,
                1.3072036044702444, 1.8764158215614224] * 1e-28
        @test all(σH_continuum.(λ, 6e3, 1e20, 4e-5) .≈ prev)
        @test argmax(σH_continuum.(λ, 6e3, 1e20, 4e-5)) == 3
        @test σH_atoms_bf(σ_atom_itp_empty, atoms_empty, 500., 6000., 1e20) == 0.0
        # Check that cross section matches data from atom with no stimulated emission
        Ly_cont = atoms[2].continua[1]  # Lyman continuum
        @test σH_atoms_bf(σ_atom_itp, atoms, Ly_cont.λ[end], 0., 1e20) ≈ Ly_cont.σ[end]
        @test σH_atoms_bf(σ_atom_itp, atoms, λ[1], 6e3, 1e20) ≈ 2.8831654074700022e-31
        @test σH_atoms_bf(σ_atom_itp, atoms, λ[end], 6e3, 1e20) ≈ 4.787807147547797e-31
    end

    @testset "σ_atoms_bf_tables" begin
        @test length(σ_atom_itp) == length(atoms)
        @test length(σ_atom_itp_empty) == 0
        @test length.(σ_atom_itp) == length.([atom.continua for atom in atoms])
        # Test against values from atomfile * abundances:
        # Using current abundances:
        @test σ_atom_itp[1][1](50.427) / abund[atoms[1].element] == 6.783e-22
        @test σ_atom_itp[1][2](14.609) / abund[atoms[1].element] == 4.844e-23
        @test σ_atom_itp[2][2](364.70521515693) / abund[atoms[2].element] ≈ 1.3788716766e-21
        # Using previous abundances:
        @test σ_atom_itp[1][1](50.427) ≈ 6.783e-22 * 0.08203515443298176
        @test σ_atom_itp[1][2](14.609) ≈ 4.844e-23 * 0.08203515443298176
        @test σ_atom_itp[2][2](364.70521515693) ≈ 1.3788716752197066e-21
        # A few edges:
        @test σ_atom_itp[1][1](50.428) == 0
        @test σ_atom_itp[1][1](12.7) == 0
        @test σ_atom_itp[2][1](22.793) == 0
        @test σ_atom_itp[2][2](364.8) == 0
        @test σ_atom_itp[2][2](91.175) == 0
    end
    @testset "σH_atoms_bf" begin
        # Against previous implementation:
        @test σH_atoms_bf(σ_atom_itp, atoms, 251.0, 6000.0, 1.0e19) ≈ 5.261436909864e-30
        @test σH_atoms_bf(σ_atom_itp, atoms, 10.0, 6000.0, 1.0e19) == 0
        @test σH_atoms_bf(σ_atom_itp, atoms, 370.0, 6000.0, 1.0e19) ≈ 1.5363426195935e-31
    end
end
