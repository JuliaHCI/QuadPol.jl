using AstroImages
using QuadPol
using Photometry

hr4796_path = joinpath(@__FILE__, "../HR4796_zimpol_stokes_cube.fits")
hr4796_cube = AstroImage(hr4796_path)

Q = hr4796_cube[:, :, 1]
U = hr4796_cube[:, :, 2]

pa=-118


cx, cy = size(Q) ./ 2 .+ 0.5
mask_inner = EllipticalAperture(cx, cy, 15, 117, 90-pa)
mask_outer = EllipticalAperture(cx, cy, 51, 182, 90-pa)

idxs = CartesianIndices(Q)
mask = map(idx -> mask_outer[idx] - mask_inner[idx], idxs)

# nx, ny = size(Q)
# xs = range(-nx/2 + 0.5, nx/2 - 0.5)
# ys = range(-ny/2 + 0.5, ny/2 - 0.5)
# sky_angles = atand.(-xs, ys')
# disk_angles = mod.(sky_angles .+ pa, 360)
# quad_mask = @. (disk_angles <= 45) | (disk_angles > 315)
# @. mask .*= 1 - Int(quad_mask)

Q_stamp = Q .* mask
U_stamp = U .* mask

quads = QuadPol.extract_quadrants(Q_stamp, U_stamp; pa=pa)

stats = QuadPol.calculate_statistics(quads)
QuadPol.print_statistics(stats)