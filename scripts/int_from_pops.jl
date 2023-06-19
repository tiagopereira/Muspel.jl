using Muspel
using AtomicData
using HDF5
using ProgressMeter
using Base.Threads


"""
Sample script to calculate Ca II 854.2 nm disk-centre intensity from
a RH atmosphere and populations.
"""
function do_work_rh(atmos_file, aux_file, atom_file)
    # read atomic data
    ca = read_atom(atom_file)
    my_line = ca.lines[5]  # Assumes CaII_PRD.yaml, index 5 for 854.2

    # read atmosphere from RH
    atmos = read_atmos_rh(atmos_file)
    ca_pops = h5read(aux_file, "atom_CA/populations")
    n_u = ca_pops[:, :, :, 5]
    n_l = ca_pops[:, :, :, 3]

    # Continuum opacity structures
    bckgr_atoms = [
        "Al.yaml",
        "C.yaml",
        "Ca.yaml",
        "Fe.yaml",
        "H_6.yaml",
        "He.yaml",
        "KI.yaml",
        "Mg.yaml",
        "N.yaml",
        "Na.yaml",
        "NiI.yaml",
        "O.yaml",
        "S.yaml",
        "Si.yaml",
    ]
    atom_files = [joinpath(AtomicData.get_atom_dir(), a) for a in bckgr_atoms]
    σ_itp = get_σ_itp(atmos, my_line.λ0, atom_files)

    a = LinRange(1f-4, 1f1, 20000)
    v = LinRange(-5f2, 5f2, 5000)
    voigt_itp = create_voigt_itp(a, v)

    # result
    intensity = Array{Float32, 3}(undef, my_line.nλ, atmos.nx, atmos.ny)

    p = ProgressMeter.Progress(atmos.nx)

    Threads.@threads for i in 1:atmos.nx
        buf = RTBuffer(atmos.nz, my_line.nλ, Float32)  # allocate inside for local scope
        for j in 1:atmos.ny
            calc_line_1D!(my_line, buf, atmos[j, i], n_u[:, j, i], n_l[:, j, i], σ_itp, voigt_itp)
            intensity[:, j, i] = buf.intensity
        end
        ProgressMeter.next!(p)
    end
    return intensity
end


"""
Sample script to calculate Ca II 854.2 nm disk-centre intensity from
a Multi3D atmosphere and populations.
"""
function do_work_multi3d(mesh_file, atmos_file, pops_file, atom_file)
    # read atomic data
    ca_atom = read_atom(atom_file)
    my_line = ca_atom.lines[5]  # Assumes CaII_run.yaml, index 5 for 854.2

    # read atmosphere and populations from Multi3D
    atmos = read_atmos_multi3d(mesh_file, atmos_file)
    pops = read_pops_multi3d(pops_file, atmos.nx, atmos.ny, atmos.nz, ca_atom.nlevels)
    n_u = pops[:, :, :, 5]
    n_l = pops[:, :, :, 3]

    # Continuum opacity structures
    bckgr_atoms = [
        "Al.yaml",
        "C.yaml",
        "Ca.yaml",
        "Fe.yaml",
        "H_6.yaml",
        "He.yaml",
        "KI.yaml",
        "Mg.yaml",
        "N.yaml",
        "Na.yaml",
        "NiI.yaml",
        "O.yaml",
        "S.yaml",
        "Si.yaml",
    ]
    atom_files = [joinpath(AtomicData.get_atom_dir(), a) for a in bckgr_atoms]
    σ_itp = get_σ_itp(atmos, my_line.λ0, atom_files)

    a = LinRange(1f-4, 1f1, 20000)
    v = LinRange(-5f2, 5f2, 5000)
    voigt_itp = create_voigt_itp(a, v)

    # result
    intensity = Array{Float32, 3}(undef, my_line.nλ, atmos.nx, atmos.ny)

    p = ProgressMeter.Progress(atmos.nx)

    Threads.@threads for i in 1:atmos.nx
        buf = RTBuffer(atmos.nz, my_line.nλ, Float32)  # allocate inside for local scope
        for j in 1:atmos.ny
            calc_line_1D!(my_line, buf, atmos[j, i], n_u[:, j, i], n_l[:, j, i], σ_itp, voigt_itp)
            intensity[:, j, i] = buf.intensity
        end
        ProgressMeter.next!(p)
    end
    return intensity
end
