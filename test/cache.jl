println("Testing printing of EvalPoint")
print(EvalPoint([0.5; 10], [1.4]))

@testset "Initialization of eval point parameters" begin
    @test_throws(ErrorException, EvalPoint([], [4.0; -38.0]))
    @test_throws(ErrorException, EvalPoint([-2.0; 1.5], []))

    x = EvalPoint([1.5; -7.4], [2.4; -1.4; Inf])
    @test x.inputs == [1.5; -7.4]
    @test x.outputs == [2.4; -1.4; Inf]
end

@testset "Cache testing" begin
    @test_throws(ErrorException, Cache((-1, 2), 1))
    @test_throws(ErrorException, Cache((1, 0), 1))
    @test_throws(ErrorException, Cache((1, 2), 0))
    cache = Cache((3, 1), 10)
    @test length(cache.Vᵏ) == 10
    @test length(get_Vk(cache)) == 0

    x1 = EvalPoint([2.4; -5.3; 2.1], [10])

    is_in_cache = add_cache!(cache, x1)
    @test is_in_cache == true
    @test length(get_Vk(cache)) == 1

    xcandidate = [1.0; 2.0; -2.0; 6.0]
    @test isincache(cache, xcandidate) == false

    xcandidate = [1.0; 2.0; 3.0]
    @test isincache(cache, xcandidate) == false

    xcandidate = [2.4; -5.3; 2.1]
    @test isincache(cache, xcandidate) == true

    @test add_cache!(cache, x1) == false
    @test length(get_Vk(cache)) == 1

    x2 = EvalPoint([-22; 32; 15.5], [10; 5])
    @test add_cache!(cache, x2) == false

    x3 = EvalPoint([-22; 32; 15.5], [-23])
    @test add_cache!(cache, x3) == true
    @test length(get_Vk(cache)) == 2

    # robustness when adding in cache two times the same eval point
    @test add_cache!(cache, x3) == false
    @test length(get_Vk(cache)) == 2

    x4 = EvalPoint([-22; 32; 15.5], [-23])
    @test add_cache!(cache, x4) == false
    @test length(get_Vk(cache)) == 2

    # special case
    x5 = EvalPoint([-5.0; 4.0], [1.5])
    @test add_cache!(cache, x5) == false
    @test length(get_Vk(cache)) == 2

    # writing in file
    save_cache(cache, "simple-cache-example.txt")
    counter = 0
    open("simple-cache-example.txt", "r") do io
        for ln in eachline(io)
            counter += 1
        end
    end
    @test counter == 3
    rm("simple-cache-example.txt")
end

#=
@testset "Initialization of eval point parameters" begin
    m2d = GranularMesh(2, [0.75; 2.5])
    @test_throws(ErrorException, EvalPoint([], [4.0; -38.0], m2d))
    @test_throws(ErrorException, EvalPoint([-2.0; 1.5], [], m2d))
    @test_throws(ErrorException, EvalPoint([1.5; 2.0; 4.0], [1.0], m2d))

    x = EvalPoint([1.5; -7.4], [2.4; -1.4; Inf], m2d)
    @test x.inputs == [1.5; -7.4]
    @test x.outputs == [2.4; -1.4; Inf]
end

@testset "Mesh updating" begin
    m2d = GranularMesh(2, [3.5; 17.5])
    @test m2d.b_init == m2d.b_k == [0; 1]
    @test m2d.a_k == [5; 2]

    x = EvalPoint([2.0; -6.0], [-67], m2d)
    refine_frame_size!(x.mesh)
    @test x.mesh.b_init == [0; 1]
    @test x.mesh.b_k == [0; 1]
    @test x.mesh.a_k == [2; 1]

    # test deep copy
    @test m2d.b_init == [0; 1]
    @test m2d.b_k == [0; 1]
    @test m2d.a_k == [5; 2]

end

@testset "Cache testing" begin
    @test_throws(ErrorException, Cache((-1, 2), 1))
    @test_throws(ErrorException, Cache((1, 0), 1))
    @test_throws(ErrorException, Cache((1,2), 0))
    cache = Cache((3, 1), 10)
    @test length(cache.Vk) == 10
    @test length(get_Vk(cache)) == 0

    m3d = GranularMesh(3, [3.5; 4; 15.1])
    x1 = EvalPoint([2.4; -5.3; 2.1], [10], m3d)

    is_in_cache = add_cache!(cache, x1)
    @test is_in_cache == true
    @test length(get_Vk(cache)) == 1

    xcandidate = [1.0; 2.0; -2.0; 6.0]
    @test isincache(cache, xcandidate) == false

    xcandidate = [1.0; 2.0; 3.0]
    @test isincache(cache, xcandidate) == false

    xcandidate = [2.4; -5.3; 2.1]
    @test isincache(cache, xcandidate) == true

    @test add_cache!(cache, x1) == false
    @test length(get_Vk(cache)) == 1

    x2 = EvalPoint([-22; 32; 15.5], [10; 5], m3d)
    @test add_cache!(cache, x2) == false

    x3 = EvalPoint([-22; 32; 15.5], [-23], m3d)
    @test add_cache!(cache, x3) == true
    @test length(get_Vk(cache)) == 2

    # robustness when adding in cache two times the same eval point
    @test add_cache!(cache, x3) == false
    @test length(get_Vk(cache)) == 2

    x4 = EvalPoint([-22; 32; 15.5], [-23], m3d)
    @test add_cache!(cache, x4) == false
    @test length(get_Vk(cache)) == 2

    # special case
    x5 = EvalPoint([-5.0; 4.0], [1.5], GranularMesh(2, [1.3; 4.1]))
    @test add_cache!(cache, x5) == false
    @test length(get_Vk(cache)) == 2

    # writing in file
    save_cache(cache, "simple-cache-example.txt")
    counter = 0
    open("simple-cache-example.txt", "r") do io
        for ln in eachline(io)
            counter += 1
        end
    end
    @test counter == 3
    rm("simple-cache-example.txt")
end
=#
