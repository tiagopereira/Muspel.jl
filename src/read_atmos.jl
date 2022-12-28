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
Reads atmosphere in the input format of MULTI3D. Returns a dictionary.
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
    ne = Mmap.mmap(fobj, Array{FloatT, 3}, shape) ./ u_l^3
    temperature = Mmap.mmap(fobj, Array{FloatT, 3}, shape, block_size)
    vz = Mmap.mmap(fobj, Array{FloatT, 3}, shape, block_size * 4) .* u_v
    nH = Mmap.mmap(fobj, Array{FloatT, 3}, shape, block_size * 5)  # reads rho
    nH = nH ./ (grph * u_l^3)
    close(fobj)
    return Atmosphere(
        x,
        y,
        z,
        permutedims(temperature, (3, 2, 1)),
        permutedims(vz, (3, 2, 1)),
        permutedims(ne, (3, 2, 1)),
        permutedims(reshape(nH, (1, nx, ny, nz)), (4, 3, 2, 1)),
    )
end


# For testing purposes
function read_atmos_multi3d_double(par_file, atmos_file; FloatT=Float32, grph=2.380491f-24)
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
    ne = Mmap.mmap(fobj, Array{FloatT, 3}, shape) ./ u_l^3
    temperature = Mmap.mmap(fobj, Array{FloatT, 3}, shape, block_size)
    vz = Mmap.mmap(fobj, Array{FloatT, 3}, shape, block_size * 4) .* u_v
    nH = Mmap.mmap(fobj, Array{FloatT, 3}, shape, block_size * 5)  # reads rho
    nH = nH ./ (grph * u_l^3)
    close(fobj)
    return Atmosphere(
        Float64.(x),
        Float64.(y),
        Float64.(z),
        Float64.(permutedims(temperature, (3, 2, 1))),
        Float64.(permutedims(vz, (3, 2, 1))),
        Float64.(permutedims(ne, (3, 2, 1))),
        Float64.(permutedims(reshape(nH, (1, nx, ny, nz)), (4, 3, 2, 1))),
    )
end


"""
Reads NLTE populations from MULTI3D output
"""
function read_pops_multi3d(pop_file, nx, ny, nz, nlevels; FloatT=Float32)::Array{FloatT, 4}
    u_l = ustrip(1f0u"cm" |> u"m")
    fobj = open(pop_file, "r")
    pops = Mmap.mmap(fobj, Array{FloatT, 4}, (nx, ny, nz, nlevels)) ./ u_l^3   # NLTE pops
    close(fobj)
    return permutedims(pops, (3, 2, 1, 4))
end

# For testing purposes
function read_pops_multi3d_double(pop_file, nx, ny, nz, nlevels; FloatT=Float32)::Array{Float64, 4}
    u_l = ustrip(1f0u"cm" |> u"m")
    fobj = open(pop_file, "r")
    pops = Mmap.mmap(fobj, Array{FloatT, 4}, (nx, ny, nz, nlevels)) ./ u_l^3   # NLTE pops
    close(fobj)
    return Float64.(permutedims(pops, (3, 2, 1, 4)))
end
