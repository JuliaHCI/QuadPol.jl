module QuadPol
using Printf

export quadpol

function extract_quadrants(
    Q::AbstractMatrix,
    U::AbstractMatrix;
    pa,
)
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

    return (;
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


function calculate_statistics(quadpols)
    # Quadrant sums
    Q000 = sum(quadpols.Q000)
    Q090 = sum(quadpols.Q090)
    Q180 = sum(quadpols.Q180)
    Q270 = sum(quadpols.Q270)
    U045 = sum(quadpols.U045)
    U135 = sum(quadpols.U135)
    U225 = sum(quadpols.U225)
    U315 = sum(quadpols.U315)

    Q_int = sum((Q000, Q090, Q180, Q270))
    U_int = sum((U045, U135, U225, U315))
    Q_abs_int = sum(abs, (Q000, Q090, Q180, Q270))
    U_abs_int = sum(abs, (U045, U135, U225, U315))

    Qphi_int = sum(quadpols.Qphi)
    Uphi_int = sum(quadpols.Uphi)

    delta_090_270 = (abs(Q090) - abs(Q270)) / (abs(Q090) + abs(Q270))
    delta_045_315 = (abs(U045) - abs(U315)) / (abs(U045) + abs(U315))
    delta_135_225 = (abs(U135) - abs(U225)) / (abs(U135) + abs(U225))

    lambda_000_180 = abs(Q000) / abs(Q180)
    lambda_045_135 = abs(U045) / abs(U135)
    lambda_315_225 = abs(U315) / abs(U225)

    lambda_a = (abs(Q090) + abs(Q270)) / (2 * abs(Q180))
    lambda_b = (abs(Q090) + abs(Q270)) / (abs(U135) + abs(U225))

    return (;
        Q000,
        Q090,
        Q180,
        Q270,
        U045,
        U135,
        U225,
        U315,
        Q_int,
        U_int,
        Q_abs_int,
        U_abs_int,
        Qphi_int,
        Uphi_int,
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
    @printf("Relative quadrant polarization statistics\n")
    @printf("Q000/Qphi =\t%f\tU045/Qphi = \t%f\n", stats.Q000 / stats.Qphi_int, stats.U045 / stats.Qphi_int)
    @printf("Q090/Qphi =\t%f\tU135/Qphi = \t%f\n", stats.Q090 / stats.Qphi_int, stats.U135 / stats.Qphi_int)
    @printf("Q180/Qphi =\t%f\tU225/Qphi = \t%f\n", stats.Q180 / stats.Qphi_int, stats.U225 / stats.Qphi_int)
    @printf("Q270/Qphi =\t%f\tU315/Qphi = \t%f\n", stats.Q270 / stats.Qphi_int, stats.U315 / stats.Qphi_int)
    @printf("\nQuadrant sum statistics\n")
    @printf("ΣQxxx/Qphi =\t%f\tΣUxxx/Qphi =\t%f\t\n", stats.Q_int / stats.Qphi_int, stats.U_int / stats.Qphi_int)
    @printf("Σ|Qxxx|/Qphi =\t%f\tΣ|Uxxx|/Qphi =\t%f\t\n", stats.Q_abs_int / stats.Qphi_int, stats.U_abs_int / stats.Qphi_int)
    @printf("\nLeft-right asymmetry deviations\n")
    @printf("Δ090/270 =\t%f%%\tΔ045/315 =\t%f%%\t\n", stats.delta_090_270 * 100, stats.delta_045_315 * 100)
    @printf("\t\t\t\tΔ135/225 =\t%f%%\t\n", stats.delta_135_225 * 100)
    @printf("\nBack-front parameter ratios\n")
    @printf("Λ000/180 =\t%f\tΛ045/135 =\t%f\t\n", stats.lambda_000_180, stats.lambda_045_135)
    @printf("\t\t\t\tΛ315/225 =\t%f\t\n", stats.lambda_315_225)
    @printf("\nSpecial back-front parameter ratios\n")
    @printf("Λa = (|Q090|+|Q270|)/(2|Q180|) \t=\t\t%f\n", stats.lambda_a)
    @printf("Λb = (|Q090|+|Q270|)/(|U135|+|U225|) =\t\t%f\n", stats.lambda_b)
end


function quadpol(Q, U; kwargs...)
    quads = extract_quadrants(Q, U; kwargs...)
    stats = calculate_statistics(quads)
    print_statistics(stats)
    return stats
end

end # module
