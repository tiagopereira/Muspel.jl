using AtomicData
using Muspel
using Test
using Unitful

@testset "background.jl" begin
    abund = get_solar_abundances()
    He = read_atom("test_atoms/He_test.yaml")
    H = read_atom("test_atoms/H_test.yaml")
    H_empty = read_atom("test_atoms/H_test_empty.yaml")
    atoms = Vector{AtomicModel}([He, H, H_empty])
    σ_atom_tables = σ_atoms_bf_tables(atoms)
    λ = [i*50.0 for i in 2:201];
    nλ = length(λ)
    log_ne = 15.0:0.05:20.0;
    log_temp = 3.0:0.015:6.0;
    table = Tables_σ(λ, log_temp, log_ne, atoms)
    @testset "Tables_σ" begin
        @test_throws MethodError Tables_σ(λ, collect(log_temp), log_ne, atoms)
        @test_throws MethodError Tables_σ(λ, log_temp, collect(log_ne), atoms)
        @test length(table.table_nh) == nλ
        @test table.table_nh[1](501.0, 6001.0) ≈ 1.51064785e-29
    end
    @testset "α_cont" begin
        @test α_cont(table, 10, λ[10], 1e4, 2.0e19, 1.0e20, 1.0e20, 2.0e19) ≈ 4.2467969689e-9
    end
    @testset "σ_atoms_bf_tables" begin
        @test σ_atom_tables[1][1](50.427) / abund[atoms[1].element] == 6.783e-22
        @test σ_atom_tables[1][1](50.428) == 0
        @test σ_atom_tables[1][1](12.7) == 0
    end
    @testset "σ_atoms_bf" begin
        @test σ_atoms_bf(σ_atom_tables, atoms, 251.0, 6000.0, 1.0e19) ≈ 5.2185826e-30
        @test σ_atoms_bf(σ_atom_tables, atoms, 10.0, 6000.0, 1.0e19) == 0
        @test σ_atoms_bf(σ_atom_tables, atoms, 370.0, 6000.0, 1.0e19) == 0
    end
    @testset "α_atoms_bf" begin
        @test α_atoms_bf(σ_atom_tables, atoms, 251.0, 6000.0, 1.0e19, 1.0e20) ≈ 5.2185826e-10
        @test α_atoms_bf(σ_atom_tables, atoms, 10.0, 6000.0, 1.0e19, 1.0e20) == 0
        @test α_atoms_bf(σ_atom_tables, atoms, 370.0, 6000.0, 1.0e19, 1.0e20) == 0
    end
end
