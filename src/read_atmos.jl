"""
Functions for reading model atmospheres from various formats.
"""

using HDF5
using Mmap
using FortranFiles
using DelimitedFiles


"""
    _read_rh_var!(buf, dset, filetype, index, trailing=())

Low-level read of an RH dataset into a preallocated `Float32` buffer, with
on-the-fly conversion to `Float32`. If the dataset has one more dimension
than expected from `buf` and `trailing`, that dimension is interpreted as
time and only timestep `index` is read. `trailing` holds fixed indices just
before the time dimension (e.g. the level in `hydrogen_populations`).
"""
function _read_rh_var!(
    buf::Array{Float32},
    dset::HDF5.Dataset,
    filetype::HDF5.Datatype,
    index::Integer,
    trailing::Tuple{Vararg{Int}}=(),
)
    idx = (ntuple(_ -> Colon(), ndims(buf))..., trailing...)
    if ndims(dset) == length(idx) + 1
        idx = (idx..., index)
    end
    HDF5.generic_read!(buf, dset, filetype, Float32, idx...)
    return buf
end


function _get_rh_dataset(fid, name)
    haskey(fid, name) && return fid[name]
    throw(ArgumentError("Dataset $name not found in $(HDF5.filename(fid))"))
end


"""
    read_atmos_rh(atmos_file; index=1, read_fullv=false, read_B=false)

Reads RH atmosphere. Returns always in single precision.

# Keywords
- `index`: timestep to read, for files with multiple timesteps (default 1).
- `read_fullv`: also read `velocity_x` and `velocity_y`, so that the
  atmosphere carries all three velocity components (default: only `velocity_z`).
- `read_B`: also read the magnetic field datasets `B_x`, `B_y`, `B_z`.
"""
function read_atmos_rh(atmos_file; index=1, read_fullv=false, read_B=false)
    h5F32 = HDF5.datatype(Float32)
    h5open(atmos_file) do fid
        magnetic = attrs(fid)["has_B"][1] == 1 ? true : false
        dims = size(fid["hydrogen_populations"])
        nz, ny, nx, nhydr = dims[1], dims[2], dims[3], dims[4]
        nt = length(dims) == 5 ? dims[5] : 1
        if !(1 <= index <= nt)
            throw(ArgumentError("Invalid timestep $index, file has $nt timestep(s)"))
        end

        x = Vector{Float32}(undef, nx)
        y = Vector{Float32}(undef, ny)
        z = Vector{Float32}(undef, nz)
        temperature = Array{Float32}(undef, nz, ny, nx)
        electron_density = Array{Float32}(undef, nz, ny, nx)
        vz = Array{Float32}(undef, nz, ny, nx)
        _read_rh_var!(x, fid["x"], h5F32, index)
        _read_rh_var!(y, fid["y"], h5F32, index)
        _read_rh_var!(z, fid["z"], h5F32, index)
        _read_rh_var!(temperature, fid["temperature"], h5F32, index)
        _read_rh_var!(electron_density, fid["electron_density"], h5F32, index)
        _read_rh_var!(vz, fid["velocity_z"], h5F32, index)

        velocity = (z = vz,)
        if read_fullv
            vx = _read_rh_var!(similar(vz), _get_rh_dataset(fid, "velocity_x"), h5F32, index)
            vy = _read_rh_var!(similar(vz), _get_rh_dataset(fid, "velocity_y"), h5F32, index)
            velocity = (x = vx, y = vy, z = vz)
        end

        magnetic_field = nothing
        if read_B && magnetic
            bx = _read_rh_var!(similar(vz), _get_rh_dataset(fid, "B_x"), h5F32, index)
            by = _read_rh_var!(similar(vz), _get_rh_dataset(fid, "B_y"), h5F32, index)
            bz = _read_rh_var!(similar(vz), _get_rh_dataset(fid, "B_z"), h5F32, index)
            magnetic_field = (x = bx, y = by, z = bz)
        end

        # Read hydrogen populations level by level: avoids allocating the
        # full (nz, ny, nx, nhydr) array when only totals are needed.
        dset_pops = fid["hydrogen_populations"]
        proton_density = similar(vz)
        if nhydr == 1
            hydrogen1_density = similar(vz)
            _read_rh_var!(hydrogen1_density, dset_pops, h5F32, index, (1,))
            # compute proton density from Saha ionisation fraction
            Threads.@threads for i in eachindex(temperature)
                @inbounds begin
                    h1 = hydrogen1_density[i]
                    proton_density[i] = h1 *
                        h_ionfrac_saha(temperature[i], electron_density[i])
                    hydrogen1_density[i] = h1 - proton_density[i]
                end
            end
        elseif nhydr == 2
            hydrogen1_density = similar(vz)
            _read_rh_var!(hydrogen1_density, dset_pops, h5F32, index, (1,))
            _read_rh_var!(proton_density, dset_pops, h5F32, index, (2,))
        else
            hydrogen1_density = zeros(Float32, nz, ny, nx)
            buffer = similar(vz)
            for level in 1:(nhydr - 1)
                _read_rh_var!(buffer, dset_pops, h5F32, index, (level,))
                hydrogen1_density .+= buffer
            end
            _read_rh_var!(proton_density, dset_pops, h5F32, index, (nhydr,))
        end

        Atmosphere(
            nx,
            ny,
            nz,
            x,
            y,
            z,
            temperature,
            velocity,
            magnetic_field,
            electron_density,
            hydrogen1_density,
            proton_density,
        )
    end
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
    return Atmosphere(
        nx,
        ny,
        nz,
        Float32.(x),
        Float32.(y),
        Float32.(z),
        temperature,
        (z = vz,),
        nothing,
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
            throw(ErrorException("Could not find $species populations in $aux_file"))
        else
            throw(e)
        end
    end
end


"""
Reads atmosphere in the input format of MULTI3D. Only Float32 atmospheres
are supported at the moment.
"""
function read_atmos_multi3d(mesh_file, atmos_file; grph=2.380491f-24)
    # Get parameters and height scale
    u_l = ustrip(1f0u"cm" |> u"m")
    u_v = ustrip(1f0u"km" |> u"m")
    nx::Int64, ny::Int64, nz::Int64, x, y, z = read_mesh(mesh_file)
    x .*= u_l
    y .*= u_l
    z .*= u_l
    # Get atmosphere and transpose
    fobj = open(atmos_file, "r")
    tmp = Array{Float32}(undef, nx, ny, nz)
    read!(fobj, tmp)
    electron_density = permutedims(tmp, (3, 2, 1))
    read!(fobj, tmp)
    temperature = permutedims(tmp, (3, 2, 1))
    read!(fobj, tmp)
    vx = permutedims(tmp, (3, 2, 1))
    read!(fobj, tmp)
    vy = permutedims(tmp, (3, 2, 1))
    read!(fobj, tmp)
    vz = permutedims(tmp, (3, 2, 1))
    read!(fobj, tmp)
    nH = permutedims(tmp, (3, 2, 1))
    close(fobj)
    proton_density = similar(temperature)

    # unit conversion and ion frac
    rho_to_nH = 1 / (grph * u_l^3)

    Threads.@threads for i in eachindex(temperature)
        electron_density[i] = electron_density[i] / u_l^3
        ionfrac = Muspel.h_ionfrac_saha(temperature[i], electron_density[i])
        proton_density[i] = nH[i] * rho_to_nH * ionfrac
        nH[i] *= rho_to_nH * (1 - ionfrac)
        vx[i] *= u_v
        vy[i] *= u_v
        vz[i] *= u_v
    end

    return Atmosphere(
        nx,
        ny,
        nz,
        x,
        y,
        z,
        temperature,
        (x = vx, y = vy, z = vz),
        nothing,
        electron_density,
        nH,
        proton_density,
    )
end


"""
Reads atmosphere in the input format of MULTI3D, at the same time as the
hydrogen populations. Only works for a H NLTE run. Only Float32 files
are supported at the moment.
"""
function read_atmos_hpops_multi3d(
        mesh_file, atmos_file, hpops_file;
        nlevels=6, grph=2.380491f-24
)
    # Get parameters and height scale
    u_l = ustrip(1f0u"cm" |> u"m")
    u_v = ustrip(1f0u"km" |> u"m")
    nx::Int64, ny::Int64, nz::Int64, x, y, z = read_mesh(mesh_file)
    x .*= u_l
    y .*= u_l
    z .*= u_l
    # Get hydrogen populations
    h_pops = Array{Float32}(undef, nx, ny, nz, nlevels)
    read!(hpops_file, h_pops)
    Threads.@threads for i in eachindex(h_pops)
        h_pops[i] = h_pops[i] / u_l^3
    end
    h1_pops = sum(h_pops[:, :, :, 1:end-1], dims=4)[:, :, :, 1]
    # Get atmosphere and transpose
    fobj = open(atmos_file, "r")
    tmp = Array{Float32}(undef, nx, ny, nz)
    read!(fobj, tmp)
    electron_density = permutedims(tmp, (3, 2, 1))
    read!(fobj, tmp)
    temperature = permutedims(tmp, (3, 2, 1))
    read!(fobj, tmp)
    vx = permutedims(tmp, (3, 2, 1))
    read!(fobj, tmp)
    vy = permutedims(tmp, (3, 2, 1))
    read!(fobj, tmp)
    vz = permutedims(tmp, (3, 2, 1))
    close(fobj)
    proton_density = permutedims(h_pops[:, :, :, end], (3, 2, 1))
    HI_density = permutedims(h1_pops, (3, 2, 1))

    Threads.@threads for i in eachindex(temperature)
        electron_density[i] /= u_l^3
        vx[i] *= u_v
        vy[i] *= u_v
        vz[i] *= u_v
    end

    atm = Atmosphere(
        nx,
        ny,
        nz,
        x,
        y,
        z,
        temperature,
        (x = vx, y = vy, z = vz),
        nothing,
        electron_density,
        HI_density,
        proton_density,
    )
    return atm, PermutedDimsArray(h_pops, (3, 2, 1, 4))
end


"""
Reads NLTE populations from MULTI3D output. Does NOT permute dims. Only Float32
files are supported at the moment.
"""
function read_pops_multi3d(pop_file, nx, ny, nz, nlevels)::Array{Float32, 4}
    u_l = ustrip(1f0u"cm" |> u"m")
    pops = Array{Float32}(undef, nx, ny, nz, nlevels)
    read!(pop_file, pops)
    Threads.@threads for i in eachindex(pops)
        pops[i] /= u_l^3
    end
    return PermutedDimsArray(pops, (3, 2, 1, 4))
end


"""
Reads mesh file from Bifrost or MULTI3D.
"""
function read_mesh(mesh_file)
    # Read all values into a single 1D array
    tmp::Vector{Float32} = Float32.(
        [a for a in vec(permutedims(readdlm(mesh_file))) if a != ""]
    )
    inc = 1
    nx = Int64(tmp[inc])
    inc += 1
    x = tmp[inc:inc + nx - 1]
    inc += nx
    ny = Int64(tmp[inc])
    inc += 1
    y = tmp[inc:inc + ny - 1]
    inc += ny
    nz = Int64(tmp[inc])
    inc += 1
    z = tmp[inc:end]
    return (nx, ny, nz, x, y, z)
end
