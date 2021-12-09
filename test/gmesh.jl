println("testing printing of GMesh")
print(GranularMesh(2, [0.5; 6.0]))
print(GranularMesh(3, [7; 1.8; 4.5]))

@testset "Initialization of mesh parameters" begin
    @test_throws(ErrorException, GranularMesh(-1, [2.0; 3.0]))
    @test_throws(MethodError, GranularMesh(0, []))
    @test_throws(ErrorException, GranularMesh(2, [0.0; 1.0]))
    @test_throws(ErrorException, GranularMesh(3, [1.0; 20.0; -1.0]))
    @test_throws(ErrorException, GranularMesh(5, [1.0; 2.0; 4.0; 5.0; 10.0; 4.0]))

    m3d = GranularMesh(3, [2.0 ; 0.5; 90.0])
    @test m3d.b_init == m3d.bᵏ == [0; 0; 1]
    @test m3d.aᵏ == [2; 1; 5]

    m4d = GranularMesh(4, [5.0; 1500.0; 0.02; 25.0])
    @test m4d.b_init == m4d.b_init == [0; 3; -1; 1]
    @test m4d.aᵏ == [5; 2; 1; 2]
end

@testset "Frame and poll size parameters" begin
    m1d = GranularMesh(1, [4.5])
    @test m1d.b_init == [0]
    @test m1d.aᵏ == [5]
    @test get_frame_size_parameter(m1d) == [5.0]
    @test get_mesh_size_parameter(m1d) == [1.0]

    m3d = GranularMesh(3, [16.0; 0.004; 450.0])
    @test m3d.b_init == m3d.bᵏ == [1; -2; 2]
    @test m3d.aᵏ == [2; 1; 5]
    @test get_mesh_size_parameter(m3d) == [10.0; 0.01; 100.0]
    @test get_frame_size_parameter(m3d) == [20.0; 0.01; 500.0]

end

@testset "Updating meshes" begin
    m2d = GranularMesh(2, [3.5; 17.5])
    @test m2d.b_init == m2d.bᵏ == [0; 1]
    @test m2d.aᵏ == [5; 2]

    refine_frame_size!(m2d)
    @test m2d.b_init == [0; 1]
    @test m2d.bᵏ == [0; 1]
    @test m2d.aᵏ == [2; 1]

    dir = [1.0; 2.0; 5.0]
    @test_throws(ErrorException, enlarge_frame_size!(m2d, dir))

    dir = [1.0; 4.0]
    enlarge_frame_size!(m2d, dir)
    @test m2d.b_init == [0;1]
    @test m2d.bᵏ == [0;1]
    @test m2d.aᵏ == [5;2]

    enlarge_frame_size!(m2d, dir)
    @test m2d.b_init == [0;1]
    @test m2d.bᵏ == [1;1]
    @test m2d.aᵏ == [1;5]

    mlim = GranularMesh(4, 10^(-7) * ones(4,), δ_min= 10^(-6))
    @test mlim.δ_min == 10^(-6)
    @test mlim.b_init == mlim.bᵏ == [-7; -7; -7; -7]
    @test mlim.aᵏ == [1, 1, 1, 1]

    # note that as we are below the tolerance for the mesh, we return the min mesh size
    @test get_mesh_size_parameter(mlim) == [10^(-6), 10^(-6), 10^(-6), 10^(-6)]

    refine_frame_size!(mlim)
    @test mlim.b_init == mlim.bᵏ == [-7; -7; -7; -7]
    @test mlim.aᵏ == [1, 1, 1, 1]

end

@testset "Projections on the mesh" begin
    m2d = GranularMesh(2, [1.0; 1.0])

    @test_throws(ErrorException, project_on_mesh(m2d, [1.5; 3.4; 6.0], [2.4; 3.0]))
    @test_throws(ErrorException, project_on_mesh(m2d, [1.5; 2.0], [-4.5]))

    proj = project_on_mesh(m2d, [1.5; -13.4], [-1.0; 2.4])
    @test proj ≈ [1.0; -13.6]

end
