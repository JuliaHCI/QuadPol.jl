module QuadPol
using Printf
using Measurements

export quadpol

struct Quadrants{T}
    Q000::AbstractMatrix{T}
    Q090::AbstractMatrix{T}
    Q180::AbstractMatrix{T}
    Q270::AbstractMatrix{T}
    U045::AbstractMatrix{T}
    U135::AbstractMatrix{T}
    U225::AbstractMatrix{T}
    U315::AbstractMatrix{T}
    Qphi::AbstractMatrix{T}
    Uphi::AbstractMatrix{T}
end


struct QuadrantsWithErrors{T}
    Q000::AbstractMatrix{T}
    Q090::AbstractMatrix{T}
    Q180::AbstractMatrix{T}
    Q270::AbstractMatrix{T}
    U045::AbstractMatrix{T}
    U135::AbstractMatrix{T}
    U225::AbstractMatrix{T}
    U315::AbstractMatrix{T}
    Qphi::AbstractMatrix{T}
    Uphi::AbstractMatrix{T}
    Q000_err::AbstractMatrix{T}
    Q090_err::AbstractMatrix{T}
    Q180_err::AbstractMatrix{T}
    Q270_err::AbstractMatrix{T}
    U045_err::AbstractMatrix{T}
    U135_err::AbstractMatrix{T}
    U225_err::AbstractMatrix{T}
    U315_err::AbstractMatrix{T}
    Qphi_err::AbstractMatrix{T}
    Uphi_err::AbstractMatrix{T}
end

function extract_quadrants(
    Q::AbstractMatrix{T},
    U::AbstractMatrix{T};
    pa,
) where {T}
    ## Step 1
    # artifically rotate the Stokes lobes by pa
    Qrot = Q * cosd(2 * pa) - U * sind(2 * pa)
    Urot = Q * sind(2 * pa) + U * cosd(2 * pa)

    ## Step 2
    # define apertures in the disk geometry
    nx, ny = size(Q)
    xs = range(-nx/2 + 0.5, nx/2 - 0.5)
    ys = range(-ny/2 + 0.5, ny/2 - 0.5)
    sky_angles = atand.(-xs, ys')
    disk_angles = mod.(sky_angles .+ pa, 360)
    # define azimuthal stokes
    Qphi = @. -Q * cosd(2 * sky_angles) - U * sind(2 * sky_angles)
    Uphi = @. Q * sind(2 * sky_angles) - U * cosd(2 * sky_angles)

    ## Step 3
    # integrate apertures for Q
    mask = @. (disk_angles <= 45) | (disk_angles > 315)
    Q000 = mask .* Qrot

    mask = @. 45 < disk_angles <= 135
    Q090 = mask .* Qrot

    mask = @. 135 < disk_angles <= 225
    Q180 = mask .* Qrot

    mask = @. 225 < disk_angles <= 315
    Q270 = mask .* Qrot

    # integrate apertures for U
    mask = @. 0 < disk_angles <= 90
    U045 = mask .* Urot

    mask = @. 90 < disk_angles <= 180
    U135 = mask .* Urot

    mask = @. 180 < disk_angles <= 270
    U225 = mask .* Urot

    mask = @. 270 < disk_angles <= 360
    U315 = mask .* Urot

    return Quadrants(
        Q000,
        Q090,
        Q180,
        Q270,
        U045,
        U135,
        U225,
        U315,
        Qphi,
        Uphi
    )
end

function extract_quadrants(
    Q::AbstractMatrix,
    Q_err::AbstractMatrix,
    U::AbstractMatrix,
    U_err::AbstractMatrix;
    pa,
)
    ## Step 1
    # artifically rotate the Stokes lobes by pa
    Qrot = @. Q * cosd(2 * pa) - U * sind(2 * pa)
    Urot = @. Q * sind(2 * pa) + U * cosd(2 * pa)
    Qrot_err = @. hypot(Q_err * cosd(2 * pa), U_err * sind(2 * pa))
    Urot_err = @. hypot(Q_err * sind(2 * pa), U_err * cosd(2 * pa))

    ## Step 2
    # define apertures in the disk geometry
    nx, ny = size(Q)
    xs = range(-nx/2 + 0.5, nx/2 - 0.5)
    ys = range(-ny/2 + 0.5, ny/2 - 0.5)
    sky_angles = atand.(-xs, ys')
    disk_angles = @. mod(sky_angles + pa, 360)
    # define azimuthal stokes
    Qphi = @. -Q * cosd(2 * sky_angles) - U * sind(2 * sky_angles)
    Uphi = @. Q * sind(2 * sky_angles) - U * cosd(2 * sky_angles)
    Qphi_err = @. hypot(Q_err * cosd(2 * sky_angles), U_err * sind(2 * sky_angles))
    Uphi_err = @. hypot(Q_err * sind(2 * sky_angles), U_err * cosd(2 * sky_angles))

    ## Step 3
    # integrate apertures for Q
    mask = @. (disk_angles <= 45) | (disk_angles > 315)
    Q000 = mask .* Qrot
    Q000_err = mask .* Qrot_err
    
    mask = @. 45 < disk_angles <= 135
    Q090 = mask .* Qrot
    Q090_err = mask .* Qrot_err

    mask = @. 135 < disk_angles <= 225
    Q180 = mask .* Qrot
    Q180_err = mask .* Qrot_err

    mask = @. 225 < disk_angles <= 315
    Q270 = mask .* Qrot
    Q270_err = mask .* Qrot_err

    # integrate apertures for U
    mask = @. 0 < disk_angles <= 90
    U045 = mask .* Urot
    U045_err = mask .* Urot_err

    mask = @. 90 < disk_angles <= 180
    U135 = mask .* Urot
    U135_err = mask .* Urot_err

    mask = @. 180 < disk_angles <= 270
    U225 = mask .* Urot
    U225_err = mask .* Urot_err

    mask = @. 270 < disk_angles <= 360
    U315 = mask .* Urot
    U315_err = mask .* Urot_err

    return QuadrantsWithErrors(
        Q000,
        Q090,
        Q180,
        Q270,
        U045,
        U135,
        U225,
        U315,
        Qphi,
        Uphi,
        Q000_err,
        Q090_err,
        Q180_err,
        Q270_err,
        U045_err,
        U135_err,
        U225_err,
        U315_err,
        Qphi_err,
        Uphi_err
    )
end

function calculate_statistics(quadpols::Quadrants)
    # Quadrant sums
    Q000 = sum(quadpols.Q000)
    Q090 = sum(quadpols.Q090)
    Q180 = sum(quadpols.Q180)
    Q270 = sum(quadpols.Q270)
    U045 = sum(quadpols.U045)
    U135 = sum(quadpols.U135)
    U225 = sum(quadpols.U225)
    U315 = sum(quadpols.U315)

    Q_quads = (Q000, Q090, Q180, Q270)
    U_quads = (U045, U135, U225, U315)
    Qd = sum(Q_quads)
    Ud = sum(U_quads)
    Qd_abs = sum(abs, Q_quads)
    Ud_abs = sum(abs, U_quads)

    Qphi = sum(quadpols.Qphi)
    Uphi = sum(quadpols.Uphi)
    
    left_right_Qsum = abs(Q090) + abs(Q270)
    left_right_Utop = abs(U045) + abs(U315)
    left_right_Ubot = abs(U135) + abs(U225)

    delta_090_270 = (abs(Q090) - abs(Q270)) / left_right_Qsum
    delta_045_315 = (abs(U045) - abs(U315)) / left_right_Utop
    delta_135_225 = (abs(U135) - abs(U225)) / left_right_Ubot

    lambda_000_180 = abs(Q000) / abs(Q180)
    lambda_045_135 = abs(U045) / abs(U135)
    lambda_315_225 = abs(U315) / abs(U225)

    lambda_a = left_right_Qsum / (2 * abs(Q180))
    lambda_b = left_right_Qsum / left_right_Ubot

    return (;
        Q000,
        Q090,
        Q180,
        Q270,
        U045,
        U135,
        U225,
        U315,
        Qd,
        Ud,
        Qd_abs,
        Ud_abs,
        Qphi,
        Uphi,
        delta_090_270,
        delta_045_315,
        delta_135_225,
        lambda_000_180,
        lambda_045_135,
        lambda_315_225,
        lambda_a,
        lambda_b
    )
end


function calculate_statistics(quadpols::QuadrantsWithErrors)
    # Quadrant sums
    Q000 = sum(quadpols.Q000) ± sqrt(sum(quadpols.Q000_err.^2))
    Q090 = sum(quadpols.Q090) ± sqrt(sum(quadpols.Q090_err.^2))
    Q180 = sum(quadpols.Q180) ± sqrt(sum(quadpols.Q180_err.^2))
    Q270 = sum(quadpols.Q270) ± sqrt(sum(quadpols.Q270_err.^2))
    U045 = sum(quadpols.U045) ± sqrt(sum(quadpols.U045_err.^2))
    U135 = sum(quadpols.U135) ± sqrt(sum(quadpols.U135_err.^2))
    U225 = sum(quadpols.U225) ± sqrt(sum(quadpols.U225_err.^2))
    U315 = sum(quadpols.U315) ± sqrt(sum(quadpols.U315_err.^2))

    Q_quads = (Q000, Q090, Q180, Q270)
    U_quads = (U045, U135, U225, U315)
    Qd = sum(Q_quads)
    Ud = sum(U_quads)
    Qd_abs = sum(abs, Q_quads)
    Ud_abs = sum(abs, U_quads)

    Qphi = sum(quadpols.Qphi) ± sqrt(sum(quadpols.Qphi_err.^2))
    Uphi = sum(quadpols.Uphi) ± sqrt(sum(quadpols.Uphi_err.^2))
    
    left_right_Qsum = abs(Q090) + abs(Q270)
    left_right_Utop = abs(U045) + abs(U315)
    left_right_Ubot = abs(U135) + abs(U225)

    delta_090_270 = (abs(Q090) - abs(Q270)) / left_right_Qsum
    delta_045_315 = (abs(U045) - abs(U315)) / left_right_Utop
    delta_135_225 = (abs(U135) - abs(U225)) / left_right_Ubot

    lambda_000_180 = abs(Q000) / abs(Q180)
    lambda_045_135 = abs(U045) / abs(U135)
    lambda_315_225 = abs(U315) / abs(U225)

    lambda_a = left_right_Qsum / (2 * abs(Q180))
    lambda_b = left_right_Qsum / left_right_Ubot

    return (;
        Q000,
        Q090,
        Q180,
        Q270,
        U045,
        U135,
        U225,
        U315,
        Qd,
        Ud,
        Qd_abs,
        Ud_abs,
        Qphi,
        Uphi,
        delta_090_270,
        delta_045_315,
        delta_135_225,
        lambda_000_180,
        lambda_045_135,
        lambda_315_225,
        lambda_a,
        lambda_b
    )
end

function print_statistics(stats)
    println("Relative quadrant polarization statistics")
    println("Q000/Qphi =\t$(stats.Q000 / stats.Qphi)\tU045/Qphi = \t$(stats.U045 / stats.Qphi)")
    println("Q090/Qphi =\t$(stats.Q090 / stats.Qphi)\tU135/Qphi = \t$(stats.U135 / stats.Qphi)")
    println("Q180/Qphi =\t$(stats.Q180 / stats.Qphi)\tU225/Qphi = \t$(stats.U225 / stats.Qphi)")
    println("Q270/Qphi =\t$(stats.Q270 / stats.Qphi)\tU315/Qphi = \t$(stats.U315 / stats.Qphi)")
    println("\nQuadrant sum statistics")
    println("ΣQxxx/Qphi =\t$(stats.Qd / stats.Qphi)\tΣUxxx/Qphi =\t$(stats.Ud / stats.Qphi)\t")
    println("Σ|Qxxx|/Qphi =\t$(stats.Qd_abs / stats.Qphi)\tΣ|Uxxx|/Qphi =\t$(stats.Ud_abs / stats.Qphi)\t")
    println("\nLeft-right asymmetry deviations")
    println("Δ090/270 =\t$(stats.delta_090_270 * 100)%\tΔ045/315 =\t$(stats.delta_045_315 * 100)%\t")
    println("\t\t\t\tΔ135/225 =\t$(stats.delta_135_225 * 100)%\t")
    println("\nBack-front parameter ratios")
    println("Λ000/180 =\t$(stats.lambda_000_180)\tΛ045/135 =\t$(stats.lambda_045_135)\t")
    println("\t\t\t\tΛ315/225 =\t$(stats.lambda_315_225)\t")
    println("\nSpecial back-front parameter ratios")
    println("Λa = (|Q090|+|Q270|)/(2|Q180|) \t=\t\t$(stats.lambda_a)")
    println("Λb = (|Q090|+|Q270|)/(|U135|+|U225|) =\t\t$(stats.lambda_b)")
end


function quadpol(Q, U; kwargs...)
    quads = extract_quadrants(Q, U; kwargs...)
    stats = calculate_statistics(quads)
    print_statistics(stats)
    return stats
end

end # module
