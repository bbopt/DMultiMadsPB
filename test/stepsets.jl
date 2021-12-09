@testset "generate 2n poll candidates" begin
    m2d = GranularMesh(2, [1.0; 1.0])
    @test m2d.b_init == m2d.bᵏ == [0;0]
    @test m2d.aᵏ == [1;1]

    @test_throws(ErrorException, PollSet2N(m2d, [2.0; 3.0; 4.0]))
    @test_throws(ErrorException, PollSet2N(m2d, [2.0; 3.0], last_success_direction=[10.0]))

    p1 = PollSet2N(m2d, [0.0; 0.0])
    rng = MersenneTwister(1234)
    pBase1 = generate_candidates(p1, rng)
    @test pBase1[:,1:2] == -pBase1[:,3:4]

    pBase2 = generate_candidates(p1, rng)
    @test pBase2[:,1:2] == -pBase2[:,3:4]
    @test pBase1 != pBase2

    refine_frame_size!(m2d)
    rng = MersenneTwister(1234)
    pBase3 = generate_candidates(p1, rng)
    @test pBase3[:,1:2] == -pBase3[:,3:4]
    @test norm(pBase3[:,1], Inf) == norm(get_frame_size_parameter(m2d), Inf)

    rng = MersenneTwister(1234)
    pBase4 = generate_ordered_candidates(p1, rng)
    @test pBase3 == pBase4

    rng = MersenneTwister(1234)
    p2 = PollSet2N(m2d, [0.0; 0.0]; last_success_direction=[-1.0;0.0])
    pBase5 = generate_ordered_candidates(p2, rng)
    @test pBase5 != pBase4
    @test pBase5 == [-0.5 -0.02 0.02 0.5; 
                     0.02 -0.5 0.5 -0.02]

end

@testset "generate n+1 poll candidates via reduction" begin

    m2d = GranularMesh(2, [1.0; 1.0])
    @test m2d.b_init == m2d.bᵏ == [0;0]
    @test m2d.aᵏ == [1;1]

    @test_throws(ErrorException, PollSetNp1(m2d, [-0.5]))
    @test_throws(ErrorException, PollSetNp1(m2d, [-0.5; 4.0],
                                            last_success_direction=[6.0; 2.0; -8.0]))

    p1 = PollSetNp1(m2d, [1.0; -1.5])
    rng = MersenneTwister(1234)
    l1_candidates = generate_candidates(p1, rng)

    @test (-p1.poll_center .+ l1_candidates[:, 1:2]) ≈ -(-p1.poll_center .+ l1_candidates[:, 3:4])

    p2 = PollSetNp1(m2d, [-3.1; 12.2], last_success_direction=[5.0;-0.5])
    rng = MersenneTwister(1234)
    l2_candidates = generate_candidates(p2, rng)
    @test size(l2_candidates) == (2,3)
    @test l2_candidates == [-3.1 -2.1 -4.1; 
                            11.2 12.2 13.2]

    rng = MersenneTwister(1234)
    l3_candidates = generate_ordered_candidates(p2, rng)
    @test size(l3_candidates) == (2,3)
    @test l2_candidates != l3_candidates
    @test l3_candidates == [-2.1 -3.1 -4.1;
                           12.2 11.2 13.2]

end

@testset "generate 2 poll candidates" begin

    m4d = GranularMesh(4, [11.0; 8.0; 5.2; 0.01])
    @test m4d.b_init == m4d.bᵏ == [1, 0, 0, -2]
    @test m4d.aᵏ == [1, 5, 5, 1]

    @test_throws(ErrorException, PollSet2(m4d, [-0.001; 2.3]))
    @test_throws(ErrorException, PollSet2(m4d, [0.8;10.0;-3.2; 4.27],
                                          last_success_direction=[2.1; -7]))

    p = PollSet2(m4d, [10.9; -3.7;4.4;-2.3], last_success_direction=[-2.1;3.5;-1.2;7.5])
    rng = MersenneTwister(26)
    l1_candidates = generate_candidates(p, rng)
    @test size(l1_candidates) == (4,2)
    @test l1_candidates[:,1] .- p.poll_center ≈ -(l1_candidates[:,2] .- p.poll_center)

    rng = MersenneTwister(26)
    l2_candidates = generate_ordered_candidates(p, rng)
    @test size(l2_candidates) == (4,2)
    @test l2_candidates == l1_candidates

end

@testset "generate 1 poll candidate" begin

    m3d = GranularMesh(3, [14.0; 8.0; 0.03])
    @test m3d.b_init == m3d.bᵏ == [1, 0, -1]
    @test m3d.aᵏ == [1, 5, 1]

    @test_throws(ErrorException, PollSet2(m3d, [-0.001; 2.3]))
    @test_throws(ErrorException, PollSet2(m3d, [0.8;10.0;-3.2],
                                          last_success_direction=[2.1]))

    p = PollSet1(m3d, [-1.0;0.0;1.6]; last_success_direction=[1.3;-4.5;1.9])
    rng = MersenneTwister(1234)
    l1_candidates = generate_candidates(p, rng)
    @test size(l1_candidates) == (3,1)

    rng = MersenneTwister(1234)
    l2_candidates = generate_ordered_candidates(p, rng)
    @test l2_candidates == l1_candidates

end

@testset "generate speculative search candidate" begin
    m2d = GranularMesh(2, [2.0; 1.0])
    @test m2d.b_init == m2d.bᵏ == [0, 0]
    @test m2d.aᵏ == [2, 1]

    @test_throws(ErrorException, SpeculativeSearchSet(m2d, [-1000.0]))
    @test_throws(ErrorException, SpeculativeSearchSet(m2d, [10.0; -2.0],
                                                      last_success_direction=[4.3;2.1;10.0]))

    rng = MersenneTwister(1234)
    s1 = SpeculativeSearchSet(m2d, [1.5; -3.1])
    @test isempty(generate_candidates(s1, rng)) == true

    s2 = SpeculativeSearchSet(m2d, [1.5; -3.1], last_success_direction=[-0.5; 0.01])
    @test generate_candidates(s2, rng) ≈ [-0.5; -3.06]

end
