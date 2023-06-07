"""
Functions for reading model atmospheres from various formats.
"""

using HDF5
using Mmap
using FortranFiles
using DelimitedFiles


"""
Reads RH atmosphere. Returns always in single precision.
"""
function read_atmos_rh(atmos_file; index=1)
    temperature = h5read(atmos_file, "temperature", (:, :, :, index))
    electron_density = h5read(atmos_file, "electron_density", (:, :, :, index))
    hydrogen_density = h5read(atmos_file, "hydrogen_populations", (:, :, :, :, index))
    vz = h5read(atmos_file, "velocity_z", (:, :, :, index))
    z = h5read(atmos_file, "z", (:, index))
    x = h5read(atmos_file, "x")
    y = h5read(atmos_file, "y")

    nz, ny, nx, nhydr = size(hydrogen_density)
    # must define proton_density
    if nhydr == 1
        hydrogen1_density = hydrogen_density[:, :, :, 1]
        Threads.@threads for i in eachindex(temperature)
            ionfrac = h_ionfrac_saha(temperature[i], electron_density[i])
            proton_density[i] = hydrogen1_density[i] * ionfrac
            hydrogen1_density[i] *= (1 - ionfrac)
        end
    elseif nhydr == 2
        hydrogen1_density = hydrogen_density[:, :, :, 1]
        proton_density = hydrogen_density[:, :, :, 2]
    elseif nhydr > 2
        proton_density = hydrogen_density[:, :, :, end]
        hydrogen1_density = dropdims(
            sum(view(hydrogen_density, :, :, :, 1:nhydr-1), dims=4);
            dims=4
        )
    end
    return Atmosphere1D(
        nx,
        ny,
        nz,
        z,
        temperature,
        vz,
        Float32.(electron_density),
        hydrogen1_density,
        proton_density
    )
end


"""
Reads RH atmosphere. Returns always in single precision.
"""
function read_atmos_hpops_rh(atmos_file, aux_file; index=1)
    temperature = h5read(atmos_file, "temperature", (:, :, :, index))
    electron_density = h5read(atmos_file, "electron_density", (:, :, :, index))
    #hydrogen_density = h5read(atmos_file, "hydrogen_populations", (:, :, :, :, index))
    vz = h5read(atmos_file, "velocity_z", (:, :, :, index))
    z = h5read(atmos_file, "z", (:, index))
    x = h5read(atmos_file, "x")
    y = h5read(atmos_file, "y")
    hydrogen_density = read_pops_rh(aux_file, "H")
    nz, ny, nx, nhydr = size(hydrogen_density)
    proton_density = hydrogen_density[:, :, :, end]
    hydrogen1_density = dropdims(
        sum(view(hydrogen_density, :, :, :, 1:nhydr-1), dims=4);
        dims=4
    )
    return Atmosphere1D(
        nx,
        ny,
        nz,
        z,
        temperature,
        vz,
        Float32.(electron_density),
        hydrogen1_density,
        proton_density
    ), hydrogen_density
end


"""
Reads array with populations for a given species.
"""
function read_pops_rh(aux_file, species)::Array{Float32, 4}
    atom = uppercase(species)
    try
        populations = h5read(aux_file, "atom_$atom/populations")
        return populations
    catch e
        if isa(e, KeyError)
            error("Could not find $species populations in $aux_file")
        else
            throw(e)
        end
    end
end


"""
Reads atmosphere in the input format of MULTI3D.
"""
function read_atmos_multi3d(mesh_file, atmos_file; FloatT=Float32, grph=2.380491f-24)
    # Get parameters and height scale
    u_l = ustrip(1f0u"cm" |> u"m")
    u_v = ustrip(1f0u"km" |> u"m")
    nx, ny, nz, x, y, z = read_mesh(mesh_file; FloatT=FloatT)
    x .*= u_l
    y .*= u_l
    z .*= u_l
    # Get atmosphere
    fobj = open(atmos_file, "r")
    shape = (nx, ny, nz)
    temperature = Array{FloatT}(undef, nx, ny, nz)
    electron_density = similar(temperature)
    vx = similar(temperature)
    vy = similar(temperature)
    vz = similar(temperature)
    nH = similar(temperature)
    read!(fobj, electron_density)
    read!(fobj, temperature)
    read!(fobj, vx)
    read!(fobj, vy)
    read!(fobj, vz)
    read!(fobj, nH)  # rho
    close(fobj)
    rho_to_nH = 1 / (grph * u_l^3)
    # Transposed arrays
    temperature_tr = Array{FloatT}(undef, nz, nx, ny)
    electron_density_tr = similar(temperature_tr)
    proton_density_tr = similar(temperature_tr)
    vx_tr = similar(temperature_tr)
    vy_tr = similar(temperature_tr)
    vz_tr = similar(temperature_tr)
    nH_tr = similar(temperature_tr)
    # Convert units, get ionisation, transpose
    Threads.@threads for i in 1:nz
        for j in 1:ny, k in 1:nx
            ne = electron_density[k, j, i] / u_l^3
            ionfrac = Muspel.h_ionfrac_saha(temperature[k, j, i], ne)
            proton_density_tr[i, j, k] = nH[k, j, i] * rho_to_nH * ionfrac
            nH_tr[i, j, k] = nH[k, j, i] * rho_to_nH * (1 - ionfrac)
            temperature_tr[i, j, k] = temperature[k, j, i]
            electron_density_tr[i, j, k] = ne
            vx_tr[i, j, k] = vx[k, j, i] * u_v
            vy_tr[i, j, k] = vy[k, j, i] * u_v
            vz_tr[i, j, k] = vz[k, j, i] * u_v
        end
    end
    return Atmosphere3D(
        Int64(nx),
        Int64(ny),
        Int64(nz),
        x,
        y,
        z,
        temperature_tr,
        vx_tr,
        vy_tr,
        vz_tr,
        electron_density_tr,
        nH_tr,
        proton_density_tr,
    )
end


"""
Reads atmosphere in the input format of MULTI3D, at the same time as the
hydrogen populations. Only works for a H NLTE run.
"""
function read_atmos_hpops_multi3d(
        mesh_file, atmos_file, hpops_file;
        nlevels=6, FloatT=Float32, grph=2.380491f-24
)
    # Get parameters and height scale
    u_l = ustrip(1f0u"cm" |> u"m")
    u_v = ustrip(1f0u"km" |> u"m")
    nx, ny, nz, x, y, z = read_mesh(mesh_file; FloatT=FloatT)
    x .*= u_l
    y .*= u_l
    z .*= u_l
    # Get hydrogen populations
    h_pops = Array{FloatT}(undef, nx, ny, nz, nlevels)
    read!(hpops_file, h_pops)
    Threads.@threads for i in eachindex(h_pops)
        h_pops[i] = h_pops[i] / u_l^3
    end
    h1_pops = sum(h_pops[:, :, :, 1:end-1], dims=4)[:, :, :, 1]
    # Get atmosphere
    fobj = open(atmos_file, "r")
    shape = (nx, ny, nz)
    block_size = nx * ny * nz * sizeof(FloatT)
    temperature = Array{FloatT}(undef, nx, ny, nz)
    electron_density = similar(temperature)
    vx = similar(temperature)
    vy = similar(temperature)
    vz = similar(temperature)
    read!(fobj, electron_density)
    read!(fobj, temperature)
    read!(fobj, vx)
    read!(fobj, vy)
    #seek(fobj, block_size * 4)
    read!(fobj, vz)
    close(fobj)
    # Transposed arrays
    temperature_tr = Array{FloatT}(undef, nz, nx, ny)
    electron_density_tr = similar(temperature_tr)
    proton_density_tr = similar(temperature_tr)
    h1_pops_tr = similar(temperature_tr)
    vx_tr = similar(temperature_tr)
    vy_tr = similar(temperature_tr)
    vz_tr = similar(temperature_tr)
    # Convert units, transpose
    Threads.@threads for i in 1:nz
        for j in 1:ny, k in 1:nx
            vx_tr[i, j, k] = vx[k, j, i] * u_v
            vy_tr[i, j, k] = vy[k, j, i] * u_v
            vz_tr[i, j, k] = vz[k, j, i] * u_v
            h1_pops_tr[i, j, k] = h1_pops[k, j, i]
            temperature_tr[i, j, k] = temperature[k, j, i]
            proton_density_tr[i, j, k] = h_pops[k, j, i, end]
            electron_density_tr[i, j, k] = electron_density[k, j, i] / u_l^3
        end
    end
    atm = Atmosphere3D(
        Int64(nx),
        Int64(ny),
        Int64(nz),
        x,
        y,
        z,
        temperature_tr,
        vx_tr,
        vy_tr,
        vz_tr,
        electron_density_tr,
        h1_pops_tr,
        proton_density_tr,
    )
    return atm, PermutedDimsArray(h_pops, (3, 2, 1, 4))
end


"""
Reads NLTE populations from MULTI3D output. Does NOT permute dims.
"""
function read_pops_multi3d(pop_file, nx, ny, nz, nlevels; FloatT=Float32)::Array{FloatT, 4}
    u_l = ustrip(1f0u"cm" |> u"m")
    pops = Array{FloatT}(undef, nx, ny, nz, nlevels)
    read!(pop_file, pops)
    Threads.@threads for i in eachindex(pops)
        pops[i] /= u_l^3
    end
    return PermutedDimsArray(pops, (3, 2, 1, 4))
end


"""
Reads mesh file from Bifrost or MULTI3D.
"""
function read_mesh(mesh_file; FloatT=Float32)
    # Read all values into a single 1D array
    tmp = [a for a in vec(permutedims(readdlm(mesh_file))) if a != ""]
    inc = 1
    nx = Int32(tmp[inc])
    inc += 1
    x = FloatT.(tmp[inc:inc + nx - 1])
    inc += nx
    ny = Int32(tmp[inc])
    inc += 1
    y = FloatT.(tmp[inc:inc + ny - 1])
    inc += ny
    nz = Int32(tmp[inc])
    inc += 1
    z = FloatT.(tmp[inc:end])
    return (nx, ny, nz, x, y, z)
end
