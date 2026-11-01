using AtomicData

# Larmor shift in nm, with λ in nm and B in T
larmor_test(λ, B) = 4.6686e-8 * λ^2 * B

@testset "intensity.jl" begin
    FALC_RH_file = joinpath(@__DIR__, "..", "data", "atmospheres", "FALC.hdf5")
    H_ATOM_file = "test_atoms/H_test.yaml"
    atm = read_atmos_rh(FALC_RH_file)[:, 1, 1]
    @testset "calc_line_1D!" begin
        h_atom = read_atom(H_ATOM_file)
        my_line = h_atom.lines[5]  # Halpha for example
        # Populations from an RH run with the FALC atmosphere
        n_l = [1.16279f5, 1.35481f5, 1.61520f5, 2.12806f5, 3.03263f5, 3.80475f5,
               5.04064f5, 7.26205f5, 1.18963f6, 1.64719f6, 2.40792f6, 3.87856f6,
               7.30773f6, 1.61016f7, 4.53617f7, 1.80606f8, 1.06021f9, 3.95242f9,
               8.26403f9, 1.4043f10, 1.8359f10, 2.0588f10, 2.2119f10, 2.2901f10,
               2.3004f10, 2.2431f10, 2.1327f10, 2.0180f10, 1.8993f10, 1.8701f10,
               1.9338f10, 2.0760f10, 2.1978f10, 2.2854f10, 2.4134f10, 2.6021f10,
               2.9463f10, 3.5499f10, 4.5300f10, 5.8483f10, 8.1597f10, 1.1706f11,
               1.5797f11, 2.0822f11, 2.6417f11, 2.9828f11, 2.9801f11, 2.7552f11,
               2.4179f11, 2.0283f11, 1.6012f11, 1.0986f11, 6.8648f10, 4.9831f10,
               4.3881f10, 5.1496f10, 8.0354f10, 1.8207f11, 4.3557f11, 1.0726f12,
               2.6948f12, 6.8632f12, 1.1562f13, 2.0921f13, 4.2045f13, 8.8683f13,
               2.0155f14, 5.0060f14, 1.0288f15, 2.0863f15, 3.5157f15, 6.0941f15,
               1.0787f16, 2.1219f16, 4.3268f16, 8.4872f16, 1.5797f17, 2.8450f17,
               4.8905f17, 8.0779f17, 1.2209f18, 1.7534f18]
        n_u = [5.90943f4, 6.76884f4, 7.90113f4, 1.00357f5, 1.35520f5, 1.63550f5,
               2.05444f5, 2.73754f5, 3.96664f5, 5.00951f5, 6.50285f5, 8.87271f5,
               1.30742f6, 2.05471f6, 3.73978f6, 8.22859f6, 2.24261f7, 5.01449f7,
               8.01383f7, 1.13485f8, 1.35617f8, 1.46254f8, 1.53010f8, 1.55528f8,
               1.53999f8, 1.48802f8, 1.40722f8, 1.32848f8, 1.24931f8, 1.24270f8,
               1.30920f8, 1.45870f8, 1.61694f8, 1.75948f8, 1.94269f8, 2.17835f8,
               2.59706f8, 3.35659f8, 4.67752f8, 6.60090f8, 1.03705f9, 1.70055f9,
               2.60173f9, 3.90411f9, 5.56120f9, 6.93626f9, 7.36287f9, 7.09179f9,
               6.39042f9, 5.46844f9, 4.37954f9, 3.03579f9, 1.90705f9, 1.38767f9,
               1.22380f9, 1.43808f9, 2.24854f9, 5.12172f9, 1.2405f10, 3.1392f10,
               8.3078f10, 2.2850f11, 4.0413f11, 7.7443f11, 1.6682f12, 3.8039f12,
               9.4887f12, 2.6459f13, 6.0266f13, 1.3630f14, 2.5012f14, 4.7587f14,
               9.3023f14, 2.0650f15, 4.7931f15, 1.0635f16, 2.2174f16, 4.4491f16,
               8.4517f16, 1.5316f17, 2.4973f17, 3.8322f17]
        σ_itp = get_σ_itp(atm, my_line.λ0, empty([""]))  # no background atoms
        a = LinRange(1f-4, 1f1, 20000)
        v = LinRange(0f2, 5f2, 2500)
        voigt_itp = create_voigt_itp(a, v)
        buf = RTBuffer(atm.nz, my_line.nλ, Float32)
        calc_line_prep!(my_line, buf, atm, σ_itp)
        calc_line_1D!(my_line, buf, my_line.λ, atm, n_u, n_l, voigt_itp;
                      to_end=false, initial_condition=:source, calc_τ_one=true)

        # Consistency tests
        @test buf.intensity[end] == buf.int_tmp[1]
        S_Planck = blackbody_λ.(my_line.λ0, atm.temperature)
        @test isapprox(buf.j_c ./ buf.α_c, buf.source_function, rtol=1e-2)  # S = S_c
        @test all(buf.source_function .> 0)
        @test all(buf.α_c .> 0)
        @test all(buf.j_c .> 0)
        @test all(isfinite.(buf.source_function))
        @test all(isfinite.(buf.intensity))
        @test all(isfinite.(buf.int_tmp))
        @test all(isfinite.(buf.ΔλD))
        @test all(isfinite.(buf.γ))
        @test all(buf.ΔλD .< 0.1)
        @test all(buf.γ .>= my_line.γ.natural)
        # Tests for this specific line profile
        @test argmin(buf.intensity) == my_line.nλ ÷ 2 + 1  # no velocity, min at line centre
        @test buf.intensity[1] > buf.intensity[my_line.nλ ÷ 2 + 1]  # absorption line
        @test isapprox(buf.intensity[my_line.nλ ÷ 2 + 1], 4.6703334, rtol=1e-3)  # line core
        @test isapprox(buf.intensity[1], 27.27908, rtol=1e-3)   # continuum intensity
        @test isapprox(buf.intensity[19], 15.57203, rtol=1e-3)  # line wing
        @test isapprox(  # symmetric line profile
            buf.intensity[1:my_line.nλ ÷ 2],
            reverse(buf.intensity[my_line.nλ ÷ 2+2:end]),
            rtol=1e-3
        )
        @test all(diff(buf.intensity[1:my_line.nλ ÷ 2]) .< 0)        # smooth decrease
        @test all(diff(buf.intensity[my_line.nλ ÷ 2 + 2:end]) .> 0)  # smooth increase

        # now calculate extreme cases
        calc_line_1D!(my_line, buf, my_line.λ, atm, n_u, n_u, voigt_itp) # emission line
        @test buf.intensity[1] < buf.intensity[my_line.nλ ÷ 2 + 1]

        calc_line_1D!(my_line, buf, my_line.λ, atm, n_u * 0, n_u * 0, voigt_itp)  # no line
        @test all(buf.intensity .== buf.intensity[1])
        @test isapprox(buf.intensity[1], 28.538631, rtol=1e-3)

        calc_line_1D!(my_line, buf ,my_line.λ, atm, n_u, n_l, voigt_itp;
                      to_end=true, initial_condition=:zero)
        @test buf.intensity[end] == buf.int_tmp[end]

        # Test isotopes
        iso_fraction = [0.5, 0.5]
        iso_weight = [my_line.mass, my_line.mass]
        iso_Δλ = [-0.1, 0.1]
        calc_line_1D_isotopes!(my_line, buf, my_line.λ, atm, n_u, n_l, voigt_itp,
                              iso_fraction, iso_weight, iso_Δλ)
        @test buf.intensity[85] ≈ buf.intensity[argmin(buf.intensity)]
        @test buf.intensity[99÷2] > minimum(buf.intensity)
    end
    @testset "calc_line_stokes_1D!" begin
        ca_atom = read_atom("test_atoms/CaII_PRD.yaml")
        line = ca_atom.lines[end]  # Ca II 854.4 nm
        nz = atm.nz
        # LTE populations
        nH = atm.hydrogen1_density .+ atm.proton_density
        pops = [saha_boltzmann(ca_atom, atm.temperature[i], atm.electron_density[i],
                               10^(6.34 - 12) * nH[i]) for i in 1:nz]
        n_up = Float32[p[5] for p in pops]
        n_lo = Float32[p[3] for p in pops]
        σ_itp = get_σ_itp(atm, line.λ0, empty([""]))
        voigt_itp = create_voigt_itp(LinRange(1f-4, 1f1, 300), LinRange(0f0, 5f2, 5001))
        # Fine table for derivatives, only valid close to the line core
        voigt_fine = create_voigt_itp(LinRange(1f-4, 1f1, 300), LinRange(0f0, 50f0, 10001))
        buf = RTBuffer(nz, line.nλ, Float32)
        calc_line_prep!(line, buf, atm, σ_itp)
        calc_line_1D!(line, buf, line.λ, atm, n_up, n_lo, voigt_itp)
        buf_s = RTBufferStokes(nz, line.nλ)
        zero_B = zeros(Float32, nz)
        with_B(bx, by, bz) = Atmosphere(
            atm.nx, atm.ny, atm.nz, atm.x, atm.y, atm.z, atm.temperature, atm.velocity,
            (x=bx, y=by, z=bz), atm.electron_density, atm.hydrogen1_density, atm.proton_density
        )
        function stokes_for(atm_B, λ=line.λ, itp=voigt_itp)
            buf_l = RTBufferStokes(nz, length(λ))
            calc_line_prep!(line, buf_l, atm_B, σ_itp)
            calc_line_stokes_1D!(line, buf_l, λ, atm_B, n_up, n_lo, itp)
            return buf_l
        end
        B0 = 0.1f0  # T
        # No field: same as scalar solution and unpolarised
        b0 = stokes_for(with_B(zero_B, zero_B, zero_B))
        @test all(b0.ΔλB .== 0)
        @test b0.α_c == buf.α_c
        @test b0.stokes[:, 1] ≈ buf.intensity rtol=1e-4
        @test all(abs.(b0.stokes[:, 2:4]) .< 1e-6 * maximum(b0.stokes[:, 1]))
        # Longitudinal field: V only, antisymmetric with field direction, I unchanged by sign
        bl = stokes_for(with_B(zero_B, zero_B, fill(B0, nz)))
        bl_neg = stokes_for(with_B(zero_B, zero_B, fill(-B0, nz)))
        @test all(bl.cosγ .== 1) && all(bl.sin²γ .== 0)
        @test bl.ΔλB[1] ≈ larmor_test(line.λ0, B0) rtol=1e-3
        @test all(bl.stokes[:, 2:3] .== 0)
        @test maximum(abs.(bl.stokes[:, 4])) > 1e-2 * maximum(bl.stokes[:, 1])
        @test bl_neg.stokes[:, 4] ≈ -bl.stokes[:, 4]
        @test bl_neg.stokes[:, 1] ≈ bl.stokes[:, 1]
        # Transverse field: Q and U only, V = 0. Rotating the field by 90° flips Q,
        # and by 45° turns Q into U
        bx = stokes_for(with_B(fill(B0, nz), zero_B, zero_B))
        by = stokes_for(with_B(zero_B, fill(B0, nz), zero_B))
        bd = stokes_for(with_B(fill(B0 / √2f0, nz), fill(B0 / √2f0, nz), zero_B))
        Qmax = maximum(abs.(bx.stokes[:, 2]))
        @test Qmax > 1e-2 * maximum(bx.stokes[:, 1])
        @test all(abs.(bx.stokes[:, 3]) .< 1e-6 * Qmax)
        @test all(abs.(bx.stokes[:, 4]) .< 1e-6 * Qmax)
        @test by.stokes[:, 2] ≈ -bx.stokes[:, 2]
        @test all(abs.(by.stokes[:, 3]) .< 1e-5 * Qmax)
        @test all(abs.(bd.stokes[:, 2]) .< 1e-5 * Qmax)
        @test bd.stokes[:, 3] ≈ bx.stokes[:, 2] rtol=1e-4 atol=1e-5 * Qmax
        # Constant oblique field: all Stokes parameters non-zero, and I,Q,U invariant to
        # rotation around the line of sight combined with the azimuth
        bg = stokes_for(with_B(fill(0.05f0, nz), fill(0.03f0, nz), fill(0.06f0, nz)))
        @test all(isfinite.(bg.stokes))
        @test all(maximum(abs.(bg.stokes[:, k])) > 1e-3 for k in 2:4)
        @test bg.stokes[:, 1] ≈ buf.intensity rtol=0.1  # I is only weakly affected
        # Weak field approximation, V = -g_eff ΔλB dI/dλ, on a fine wavelength grid
        λf = collect(range(line.λ0 - 0.08, line.λ0 + 0.08, length=801))
        bw = stokes_for(with_B(zero_B, zero_B, fill(1f-3, nz)), λf, voigt_fine)
        I = bw.stokes[:, 1]
        V_weak = -bw.ΔλB[1] .* (I[3:end] - I[1:end-2]) ./ (λf[3:end] - λf[1:end-2])
        sel = abs.(V_weak) .> 0.1 * maximum(abs.(V_weak))
        ratio = sort(bw.stokes[2:end-1, 4][sel] ./ V_weak[sel])
        geff = sum(line.σb_strength .* line.σb_shift)
        @test geff ≈ 1.1 rtol=1e-3
        @test isapprox(ratio[end ÷ 2 + 1], geff, rtol=0.1)
        # Line without Zeeman components
        h_line = read_atom(H_ATOM_file).lines[5]
        @test_throws ArgumentError calc_line_stokes_1D!(
            h_line, buf_s, h_line.λ, with_B(zero_B, zero_B, zero_B), n_up, n_lo, voigt_itp
        )
    end
    @testset "calc_τ_cont!" begin
        σ_itp = get_σ_itp(atm, 500f0, empty([""]))
        tau1 = Vector{Float32}(undef, atm.nz)
        tau2 = Vector{Float32}(undef, atm.nz)
        calc_τ_cont!(atm, tau1, σ_itp)
        @test tau1[1] == 0
        # Testing against implementation
        @test tau1[end] ≈ 19.67434
        @test minimum(log10.(tau1[2:end])) > -9
        # Testing with reversed atmosphere, final integration should be the same
        rev_atm = deepcopy(atm)
        rev_atm.temperature[:] = reverse(atm.temperature)
        rev_atm.electron_density[:] = reverse(atm.electron_density)
        rev_atm.proton_density[:] = reverse(atm.proton_density)
        rev_atm.hydrogen1_density[:] = reverse(atm.hydrogen1_density)
        rev_atm.z[:] = reverse(atm.z)
        calc_τ_cont!(rev_atm, tau2, σ_itp)
        @test tau1[end] ≈ tau2[end]
    end
    @testset "calc_τ_one_height" begin
        nz = 20
        aa = ones(nz) * 0.1
        zz = 1.0:nz
        @test calc_τ_one_height(zz, aa) ≈ 11
        # Edge cases
        aa *= 0
        @test calc_τ_one_height(zz, aa) == nz
        aa .+= 1
        @test calc_τ_one_height(zz, aa) == 0
        @test calc_τ_one_height(1f0:nz, ones(Float32, nz)) isa Float32
        @test calc_τ_one_height(1.0:nz, ones(Float64, nz)) isa Float64
    end
end
