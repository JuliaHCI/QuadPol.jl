using AstroImages
using Photometry

@testset "HR4796" begin
    hr4796_path = joinpath(dirname(@__FILE__), "HR4796_zimpol_stokes_cube.fits")
    hr4796_cube = AstroImage(hr4796_path)
    Q = hr4796_cube[:, :, 1]
    U = hr4796_cube[:, :, 2]

    pa=-118

    cx, cy = size(Q) ./ 2 .+ 0.5
    mask_inner = EllipticalAperture(cx, cy, 15, 117, 90-pa)
    mask_outer = EllipticalAperture(cx, cy, 51, 182, 90-pa)

    idxs = map(intersect, axes(mask_outer), axes(Q)) |> CartesianIndices
    mask = map(idx -> mask_outer[idx] - mask_inner[idx], idxs)

    Q_stamp = Q[idxs] .* mask
    U_stamp = U[idxs] .* mask
    quads = QuadPol.extract_quadrants(Q_stamp, U_stamp; pa=pa)
    @testset "without errors" begin
        stats = QuadPol.calculate_statistics(quads)
        QuadPol.print_statistics(stats)
        @test stats.Qd < stats.Qd_abs
        @test stats.Ud < stats.Ud_abs
    end
    @testset "with errors" begin
        quads2 = QuadPol.extract_quadrants(Q_stamp, quads.Uphi, U_stamp, quads.Uphi; pa=pa)
        stats2 = QuadPol.calculate_statistics(quads2)
        QuadPol.print_statistics(stats2)
        @test stats2.Qd < stats2.Qd_abs
        @test stats2.Ud < stats2.Ud_abs
    end
end