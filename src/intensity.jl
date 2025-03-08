# Functions to calculate intensity and related quantities
const HC_K = ustrip((h * c_0 / k_B) |> u"K * nm")
const αCONST = ustrip((e^2 / (4*ε_0*m_e*c_0^2)) |> u"nm") / sqrt(π)
const INV_M3_to_NM3 = Float32(ustrip(1u"m^-3" |> u"nm^-3"))


"""
    function calc_line_1D!(
        line::AtomicLine,
        buf::RTBuffer{T},
        atm::Atmosphere1D{T},
        n_up::AbstractVector{T},
        n_lo::AbstractVector{T},
        σ_itp::ExtinctionItpNLTE{<:Real},
        voigt_itp::Interpolations.AbstractInterpolation{<:Number, 2};
        to_end::Bool=false,
        initial_condition=:source
    )

Calculate emerging disk-centre intensity for a given line in a 1D atmosphere.
"""
function calc_line_1D!(
    line::AtomicLine,
    buf::RTBuffer{T},
    atm::Atmosphere1D{1, T},
    n_up::AbstractVector{T},
    n_lo::AbstractVector{T},
    σ_itp::ExtinctionItpNLTE{<:Real},
    voigt_itp::Interpolations.AbstractInterpolation{<:Number, 2};
    to_end::Bool=false,
    initial_condition=:source,
    calc_τ_one::Bool=false,
) where T <: AbstractFloat
    γ_energy = ustrip((h * c_0 / (4 * π * line.λ0 * u"nm")) |> u"J")
    if to_end  # direction of integration
        end_point = atm.nz
        vsign = -1
    else
        end_point = 1
        vsign = 1
    end
    # wavelength-independent part (continuum + broadening + Doppler width)
    for i in 1:atm.nz
        buf.α_c[i] = α_cont(
            σ_itp,
            atm.temperature[i],
            atm.electron_density[i],
            atm.hydrogen1_density[i],
            atm.proton_density[i]
        )
        buf.j_c[i] = buf.α_c[i] * blackbody_λ(line.λ0, atm.temperature[i])
        buf.γ[i] = calc_broadening(
            line.γ,
            atm.temperature[i],
            atm.electron_density[i],
            atm.hydrogen1_density[i]
        )
        buf.ΔλD[i] = doppler_width(line.λ0, line.mass, atm.temperature[i])
    end

    # Calculate line opacity and intensity
    for (i, λ) in enumerate(line.λ)
        for iz in 1:atm.nz
            # Wavelength-dependent part
            a = damping(buf.γ[iz], λ, buf.ΔλD[iz])  # very small dependence on λ
            v = (λ - line.λ0 + line.λ0 * atm.velocity_z[iz] * vsign / ustrip(c_0)) / buf.ΔλD[iz]
            profile = real(voigt_itp(a, abs(v))) / (sqrt(π) * buf.ΔλD[iz])  # units nm^-1
            # Part that only multiplies by wavelength:
            α_tmp = γ_energy * profile
            j_tmp = α_tmp
            α_tmp *= n_lo[iz] * line.Blu - n_up[iz] * line.Bul
            j_tmp *= n_up[iz] * line.Aul
            α_tmp = α_tmp * 1f9 + buf.α_c[iz]   # convert α_tmp to m^-1
            j_tmp = j_tmp * 1f-3 + buf.j_c[iz]  # convert j_tmp to kW m^3 nm^-1
            buf.source_function[iz] = j_tmp / α_tmp
            buf.α_total[iz] = α_tmp
        end
        piecewise_1D_linear!(atm.z, buf.α_total, buf.source_function, buf.int_tmp;
                             to_end=to_end, initial_condition=initial_condition)
        buf.intensity[i] = buf.int_tmp[end_point]
        if calc_τ_one & !to_end
            buf.τ_one_height[i] = calc_τ_one_height2(atm.z, buf.α_total)
        end
    end
    return nothing
end



"""
    function calc_LTE_line_1D!(
        line::AtomicLine,
        buf::RTBuffer{T},
        atm::Atmosphere1D{T},
        nl_itp::Interpolations.AbstractInterpolation{<:Number, 2};
        σ_itp::ExtinctionItpNLTE{<:Real},
        voigt_itp::Interpolations.AbstractInterpolation{<:Number, 2};
        to_end::Bool=false,
        initial_condition=:source
    )

Calculate emerging disk-centre intensity for a given line in a 1D atmosphere.
Using simplifications for the LTE case.
"""
function calc_LTE_line_1D!(
    line::AtomicLine,
    buf::RTBufferLTE{T},
    atm::Atmosphere1D{1, T},
    abundance::Real,
    nl_itp::Interpolations.AbstractInterpolation{<:Number, 2},
    σ_itp::ExtinctionItpNLTE{<:Real},
    voigt_itp::Interpolations.AbstractInterpolation{<:Number, 2};
    to_end::Bool=false,
    initial_condition=:source,
) where T <: AbstractFloat
    if to_end  # direction of integration
        end_point = atm.nz
        vsign = -1
    else
        end_point = 1
        vsign = 1
    end
    # wavelength-independent part (continuum + broadening + Doppler width)
    for i in 1:atm.nz
        buf.α_c[i] = α_cont(
            σ_itp,
            atm.temperature[i],
            atm.electron_density[i],
            atm.hydrogen1_density[i],
            atm.proton_density[i]
        )
        buf.source_function[i] = blackbody_λ(line.λ0, atm.temperature[i])
        buf.γ[i] = calc_broadening(
            line.γ,
            atm.temperature[i],
            atm.electron_density[i],
            atm.hydrogen1_density[i]
        )
        buf.ΔλD[i] = doppler_width(line.λ0, line.mass, atm.temperature[i])
        # units: nm^-3
        n_l = nl_itp(
                atm.temperature[i],
                atm.electron_density[i]
            ) * abundance * ((atm.hydrogen1_density[i]+atm.proton_density[i]) * INV_M3_to_NM3)
        # units: nm^-1
        buf.α_l[i] = (
            αCONST * n_l * line.f_value * line.λ0^2 / buf.ΔλD[i] *
            (1 - exp(-HC_K / (atm.temperature[i]*line.λ0)))
        )
    end

    # Calculate line opacity and intensity
    for (i, λ) in enumerate(line.λ)
        for iz in 1:atm.nz
            # Wavelength-dependent part
            a = damping(buf.γ[iz], λ, buf.ΔλD[iz])  # very small dependence on λ
            v = (λ - line.λ0 + line.λ0 * atm.velocity_z[iz] * vsign / ustrip(c_0)) / buf.ΔλD[iz]
            profile = real(voigt_itp(a, abs(v)))
            # convert α_l from nm^-1 to m^-1
            buf.α_total[iz] = buf.α_c[iz] + buf.α_l[iz] * profile * 1f9
        end
        piecewise_1D_linear!(atm.z, buf.α_total, buf.source_function, buf.int_tmp;
                             to_end=to_end, initial_condition=initial_condition)
        buf.intensity[i] = buf.int_tmp[end_point]
    end
    return nothing
end


"""
Calculate continuum optical depth in the vertical direction,
from the observer to the stellar interior. The wavelength
is defined by σ_itp.
"""
function calc_τ_cont!(
    atm::Atmosphere1D{1, T},
    τ::AbstractVector{<:Real},
    σ_itp::ExtinctionItpNLTE{<:Real},
) where T <: AbstractFloat
    τ[1] = zero(T)
    α = α_cont(
        σ_itp,
        atm.temperature[1],
        atm.electron_density[1],
        atm.hydrogen1_density[1],
        atm.proton_density[1]
    )
    for i in 2:atm.nz
        α_next = α_cont(
            σ_itp,
            atm.temperature[i],
            atm.electron_density[i],
            atm.hydrogen1_density[i],
            atm.proton_density[i]
        )
        τ[i] = τ[i-1] + abs(atm.z[i] - atm.z[i-1]) * (α + α_next) / 2
        α = α_next
    end
    return nothing
end

"""
Calculate height where τ=1 using linear interpolation
"""
function calc_τ_one_height(
    height::AbstractVector{<:Real},
    α::AbstractVector{T},
) where T <: Real
    nz = length(height)
    @assert length(α) == nz "Height and extinction have different sizes"
    τ_prev = zero(T)
    τ_curr = zero(T)
    τ_one = zero(T)
    for i in 2:nz
        # Manually integrate optical depth
        tmp = abs(height[i-1] - height[i])
        τ_prev = τ_curr
        τ_curr = τ_prev + tmp * (α[i] + α[i-1]) / 2
        # Limit cases, reached at top or not at all
        if ((τ_curr > 1) & (i == 2)) | ((τ_curr < 1) & (i == nz))
            τ_one = height[i]
        elseif (τ_curr > 1) & (τ_prev < 1)
            # Manual interpolation
            τ_one = height[i-1] + (1 - τ_prev) / (τ_curr - τ_prev) * tmp
            break
        end
    end
    return τ_one
end


"""
Calculate height where τ=1. Using higher order.
"""
function calc_τ_one_height2(
    height::AbstractVector{<:Real},
    α::AbstractVector{T},
) where T <: Real
    nz = length(height)
    @assert length(α) == nz "Height and extinction have different sizes"
    τ = zeros(T, nz)
    for i in 2:nz # Manually integrate optical depth
        τ[i] = τ[i-1] + abs(height[i-1] - height[i]) * (α[i] + α[i-1]) / 2
    end
    itp = extrapolate(interpolate(τ, height, FritschCarlsonMonotonicInterpolation()), Flat())
    return itp(one(T))
end
