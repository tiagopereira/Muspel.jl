using Muspel
using BifrostTools
using AtomicData
using Base.Threads

"""
Calculates τ500 from a 3D atmosphere.
"""
function τ_from_atmos(atmos::AbstractAtmosphere{3}; wave=500)
    
    # Background atoms to source bound-free edges as sources of continuum opacity
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
    σ_itp = get_σ_itp(atmos, wave, atom_files)
    τ = similar(atmos.temperature)

    @threads for i in 1:atmos.nx
        buf = similar(atmos.z)
        for j in 1:atmos.ny
            calc_τ_cont!(atmos[:, j, i], buf, σ_itp)
            τ[:, j, i] .= buf
        end
    end
    return τ
end

"""
Calculates continuum optical depth at 500 nm for a given Bifrost snapshot. 
Assumes LTE.

Returns τ with shape (z,y,x) by default, which is the standard order for 
Muspel. `permute=true` rearranges dimensions back to (x,y,z).
"""
function τ500_from_snap(
        xp::BifrostExperiment, 
        snap::Int;
        slicex::AbstractVector{<:Integer}=Int[], 
        slicey::AbstractVector{<:Integer}=Int[],
        slicez::AbstractVector{<:Integer}=Int[],
        permute::Bool=false
)
    
    # grams per H atom
    grph = 2.38049f-24
    
    temperature = get_var(xp,snap,"tg",units="si",slicex=slicex,slicey=slicey,slicez=slicez)
    electron_density = get_electron_density(xp,snap,units="si",slicex=slicex,slicey=slicey,slicez=slicez,verbose=false)
    rho = get_var(xp,snap,"r",units="si",slicex=slicex,slicey=slicey,slicez=slicez)

    # Fix indices
    isempty(slicex) && ( slicex = 1:xp.mesh.mx )
    isempty(slicey) && ( slicey = 1:xp.mesh.my )
    isempty(slicez) && ( slicez = 1:xp.mesh.mz )

    x = xp.mesh.x[slicex] * 1f6 # to SI units
    y = -xp.mesh.y[slicey] * 1f6 # rotate and to SI units
    z = -xp.mesh.z[slicez] * 1f6 # rotate and to SI units

    nx=length(x); ny=length(y); nz=length(z)
    
    # Permute dimensions
    temperature = permutedims(temperature, (3,2,1))    
    electron_density = permutedims(electron_density, (3,2,1))
    rho = permutedims(rho, (3,2,1))
    
    hydrogen1_density = rho/(grph * 1f-3) # In units 1/m^3
    proton_density = similar(hydrogen1_density)
    @threads for i in eachindex(temperature)
        ionfrac = Muspel.h_ionfrac_saha(temperature[i], electron_density[i])
        proton_density[i] = hydrogen1_density[i] * ionfrac
        hydrogen1_density[i] *= (1 - ionfrac)
    end

    atmos = Atmosphere(
        nx,
        ny,
        nz,
        x,
        y,
        z,
        temperature,
        NamedTuple(),  # velocity not needed
        nothing,       # no magnetic field
        electron_density,
        hydrogen1_density,  # neutral hydrogen across all levels
        proton_density,
    )

    τ = τ_from_atmos(atmos, wave=500)
    if permute
        τ = permute_dims(τ, (3,2,1))
    end

    return τ
end
