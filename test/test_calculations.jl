
@testset "Trivial Quadrants" begin
    Q_img = Matrix{Float64}(undef, 50, 50)
    U_img = Matrix{Float64}(undef, 50, 50)

    nx, ny = size(Q_img)
    xs = range(-nx/2 + 0.5, nx/2 - 0.5)
    ys = range(-ny/2 + 0.5, ny/2 - 0.5)
    sky_angles = mod.(atand.(-xs, ys'), 360)

    quad_values = [1, 2, 3, 4]

    @. Q_img[(sky_angles < 45) | (315 <= sky_angles)] = quad_values[1]
    @. Q_img[(45 <= sky_angles < 135)] = quad_values[2]
    @. Q_img[(135 <= sky_angles < 225)] = quad_values[3]
    @. Q_img[(225 <= sky_angles < 315)] = quad_values[4]
    
    @. U_img[(0 <= sky_angles < 90)] = quad_values[1]
    @. U_img[(90 <= sky_angles < 180)] = quad_values[2]
    @. U_img[(180 <= sky_angles < 270)] = quad_values[3]
    @. U_img[(270 <= sky_angles < 360)] = quad_values[4]

    pa = 0
    quad = QuadPol.extract_quadrants(Q_img, U_img; pa=pa)

    @test maximum(quad.Q000) ≈ quad_values[1]
    @test maximum(quad.Q090) ≈ quad_values[2]
    @test maximum(quad.Q180) ≈ quad_values[3]
    @test maximum(quad.Q270) ≈ quad_values[4]

    @test maximum(quad.U045) ≈ quad_values[1]
    @test maximum(quad.U135) ≈ quad_values[2]
    @test maximum(quad.U225) ≈ quad_values[3]
    @test maximum(quad.U315) ≈ quad_values[4]

    stats = QuadPol.calculate_statistics(quad)
    quad_area = prod(size(Q_img)) / 4
    @test stats.Q000 ≈ quad_area * quad_values[1]
    @test stats.Q090 ≈ quad_area * quad_values[2]
    @test stats.Q180 ≈ quad_area * quad_values[3]
    @test stats.Q270 ≈ quad_area * quad_values[4]

    @test stats.U045 ≈ quad_area * quad_values[1]
    @test stats.U135 ≈ quad_area * quad_values[2]
    @test stats.U225 ≈ quad_area * quad_values[3]
    @test stats.U315 ≈ quad_area * quad_values[4]

    @test stats.Qd ≈ sum(quad_area .* quad_values)
    @test stats.Ud ≈ sum(quad_area .* quad_values)

    @test stats.Qd_abs ≈ sum(quad_area .* quad_values)
    @test stats.Ud_abs ≈ sum(quad_area .* quad_values)

    @test stats.delta_090_270 ≈ (quad_values[2] - quad_values[4]) / (quad_values[2] + quad_values[4])
    @test stats.delta_045_315 ≈ (quad_values[1] - quad_values[4]) / (quad_values[1] + quad_values[4])
    @test stats.delta_135_225 ≈ (quad_values[2] - quad_values[3]) / (quad_values[2] + quad_values[3])

    @test stats.lambda_000_180 ≈ quad_values[1]/quad_values[3]
    @test stats.lambda_045_135 ≈ quad_values[1]/quad_values[2]
    @test stats.lambda_315_225 ≈ quad_values[4]/quad_values[3]

    @test stats.lambda_a ≈ (quad_values[2] + quad_values[4]) / (2 * quad_values[3])
    @test stats.lambda_b ≈ (quad_values[2] + quad_values[4]) / (quad_values[2] + quad_values[3])    
end