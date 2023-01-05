"""
Functions for reading model atmospheres from various formats.
"""

using HDF5
using Mmap
using FortranFiles


"""
Reads RH atmosphere. Returns always in single precision.
"""
function read_atmos_rh(atmos_file; index=1)::Atmosphere{Float32}
    temperature = h5read(atmos_file, "temperature", (:, :, :, index))
    ne = h5read(atmos_file, "electron_density", (:, :, :, index))
    nH = h5read(atmos_file, "hydrogen_populations", (:, :, :, :, index))
    vz = h5read(atmos_file, "velocity_z", (:, :, :, index))
    z = h5read(atmos_file, "z", (:, index))
    x = h5read(atmos_file, "x")
    y = h5read(atmos_file, "y")
    return Atmosphere(Float32.(x), Float32.(y), z, temperature, vz, Float32.(ne), nH)
end

# For testing purposes
function read_atmos_rh_double(atmos_file; index=1)::Atmosphere{Float64}
    temperature = h5read(atmos_file, "temperature", (:, :, :, index))
    ne = h5read(atmos_file, "electron_density", (:, :, :, index))
    nH = h5read(atmos_file, "hydrogen_populations", (:, :, :, :, index))
    vz = h5read(atmos_file, "velocity_z", (:, :, :, index))
    z = h5read(atmos_file, "z", (:, index))
    x = h5read(atmos_file, "x")
    y = h5read(atmos_file, "y")
    return Atmosphere(
        Float64.(x),
        Float64.(y),
        Float64.(z),
        Float64.(temperature),
        Float64.(vz),
        Float64.(ne),
        Float64.(nH),
    )
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
This version does not permute dims. This is not currently working,
as it is not possible to reassign the hydrogen populations.
"""
function read_atmos_multi3d(par_file, atmos_file; FloatT=Float32, grph=2.380491f-24)
    # Get parameters and height scale
    u_l = ustrip(1f0u"cm" |> u"m")
    u_v = ustrip(1f0u"km" |> u"m")
    fobj = FortranFile(par_file, "r")
    _ = read(fobj, Int32)
    nx = read(fobj, Int32)
    ny = read(fobj, Int32)
    nz = read(fobj, Int32)
    x = FloatT.(read(fobj, (Float64, nx))) * u_l
    y = FloatT.(read(fobj, (Float64, ny))) * u_l
    z = FloatT.(read(fobj, (Float64, nz))) * u_l
    close(fobj)
    # Get atmosphere
    fobj = open(atmos_file, "r")
    shape = (nx, ny, nz)
    block_size = nx * ny * nz * sizeof(FloatT)
    temperature = Array{FloatT}(undef, nx, ny, nz)
    electron_density = similar(temperature)
    vz = similar(temperature)
    nH = similar(temperature)
    proton_density = similar(temperature)
    read!(fobj, electron_density)
    read!(fobj, temperature)
    # skip vx, vy
    seek(fobj, block_size * 4)
    read!(fobj, vz)
    read!(fobj, nH)  # rho
    close(fobj)
    rho_to_nH = 1 / (grph * u_l^3)
    # convert units, get hydrogen ionisation
    Threads.@threads for i in eachindex(temperature)
        electron_density[i] /= u_l^3
        vz[i] *= u_v
        ionfrac = h_ionfrac_saha(temperature[i], electron_density[i])
        proton_density[i] = nH[i] * rho_to_nH * ionfrac
        nH[i] *= rho_to_nH * (1 - ionfrac)
    end
    return AtmosphereM3D(
        Int64(nx),
        Int64(ny),
        Int64(nz),
        z,
        temperature,
        vz,
        electron_density,
        nH,
        proton_density,
    )
end


"""
Reads atmosphere in the input format of MULTI3D, at the same time as the
hydrogen populations. Only works for a H NLTE run.
This version does not permute dims.
"""
function read_atmos_hpops_multi3d(
        par_file, atmos_file, hpops_file;
        nlevels=6, FloatT=Float32, grph=2.380491f-24
)
    # Get parameters and height scale
    u_l = ustrip(1f0u"cm" |> u"m")
    u_v = ustrip(1f0u"km" |> u"m")
    fobj = FortranFile(par_file, "r")
    _ = read(fobj, Int32)
    nx = Int64(read(fobj, Int32))
    ny = Int64(read(fobj, Int32))
    nz = Int64(read(fobj, Int32))
    x = FloatT.(read(fobj, (Float64, nx))) * u_l
    y = FloatT.(read(fobj, (Float64, ny))) * u_l
    z = FloatT.(read(fobj, (Float64, nz))) * u_l
    close(fobj)
    # Get hydrogen populations
    h_pops = Array{FloatT}(undef, nx, ny, nz, nlevels)
    read!(hpops_file, h_pops)
    h_pops ./= u_l^3
    h1_pops = sum(h_pops[:, :, :, 1:end-1], dims=4)[:, :, :, 1]
    # Get atmosphere
    fobj = open(atmos_file, "r")
    shape = (nx, ny, nz)
    block_size = nx * ny * nz * sizeof(FloatT)
    temperature = Array{FloatT}(undef, nx, ny, nz)
    ne = similar(temperature)
    vz = similar(temperature)
    read!(fobj, ne)
    read!(fobj, temperature)
    # skip vx, vy
    seek(fobj, block_size * 4)
    read!(fobj, vz)
    close(fobj)
    Threads.@threads for i in eachindex(ne)
        ne[i] /= u_l^3
        vz[i] *= u_v
    end
    atm = AtmosphereM3D(nx, ny, nz, z, temperature, vz, ne, h1_pops, h_pops[:, :, :, end])
    return atm, h_pops
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
    return pops
end


# For testing purposes
function read_pops_multi3d_double(pop_file, nx, ny, nz, nlevels; FloatT=Float32)::Array{Float64, 4}
    u_l = ustrip(1f0u"cm" |> u"m")
    pops = Array{FloatT}(undef, nx, ny, nz, nlevels)
    read!(pop_file, pops)
    pops ./= u_l^3
    return Float64.(pops)
end
