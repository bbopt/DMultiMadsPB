@testset "Initialization of barrier" begin

    @test_throws(ErrorException, Barrier((0, 1), 10))
    @test_throws(ErrorException, Barrier((-3, 1), 10))
    @test_throws(ErrorException, Barrier((1, 0), 10))
    @test_throws(ErrorException, Barrier((1, -3), 10))
    @test_throws(ErrorException, Barrier((4, 3), 0))
    @test_throws(ErrorException, Barrier((5, 2), 10, h_max=-1.3))

    b1 = Barrier((5, 2), 10, h_max=0.0)
    @test b1.dims == (5, 2)
    @test b1.h_max == 0
    @test isempty(get_Fk(b1)) == true
    @test isempty(get_Uk(b1)) == true
    @test isempty(get_Ik(b1)) == true
    @test b1.max_size == length(b1.elements) == length(b1.within_Fk) == length(b1.within_Ik) == length(b1.within_Uk)
    @test b1.max_size == length(b1.meshes)

    b2 = Barrier((1, 5), 2)
    @test b2.dims == (1, 5)
    @test b2.h_max == Inf
    @test isempty(get_Fk(b2)) == true
    @test isempty(get_Uk(b2)) == true
    @test isempty(get_Ik(b2)) == true
    @test b2.max_size == length(b2.elements) == length(b2.within_Fk) == length(b2.within_Ik) == length(b2.within_Uk)
    @test b2.max_size == length(b2.meshes)
    @test b2.last_index == 0

end

@testset "Feasible insertion for two objectives" begin
    b2obj = Barrier((2, 2), 6)
    @test b2obj.max_size == length(b2obj.elements) == length(b2obj.within_Fk) == length(b2obj.within_Ik) == length(b2obj.within_Uk)
    @test b2obj.max_size == length(b2obj.meshes) == length(b2obj.parent_indexes)

    m2d = GranularMesh(2, [3.5; 17.5])
    @test_throws(ErrorException, add_feasible!(b2obj, OVector([1.0; -4.0; 6.0], 0), m2d))
    @test_throws(ErrorException, add_feasible!(b2obj, OVector([1.0; -4.0], 2.6), m2d))

    insert_flag = add_feasible!(b2obj, OVector([2.0; 3.0], 0), m2d)
    @test insert_flag == :improves
    @test length(get_Fk(b2obj)) == 1
    @test get_Fk(b2obj)[1].f == b2obj.elements[1].f ≈ [2.0; 3.0]

    # first tentative
    insert_flag = add_feasible!(b2obj, OVector([3.0; 3.0], 0), m2d)
    @test insert_flag == :is_dominated
    @test length(get_Fk(b2obj)) == 1
    @test get_Fk(b2obj)[1].f == b2obj.elements[1].f ≈ [2.0; 3.0]
    @test b2obj.elements[2].f ≈ [3.0; 3.0] && b2obj.elements[2].h == 0

    # second tentative
    insert_flag = add_feasible!(b2obj, OVector([3.0; 2.0], 0), m2d)
    @test insert_flag == :extends
    @test length(get_Fk(b2obj)) == 2
    @test get_Fk(b2obj)[2].f == b2obj.elements[3].f ≈ [3.0; 2.0]
    @test parentindices(get_Fk(b2obj))[1][2] == 3
    @test b2obj.elements[2].f ≈ [3.0; 3.0] && b2obj.elements[2].h == 0

    # third point
    insert_flag = add_feasible!(b2obj, OVector([5.0; 1.0], 0), m2d)
    @test insert_flag == :extends
    @test length(get_Fk(b2obj)) == 3
    @test get_Fk(b2obj)[3].f == b2obj.elements[4].f ≈ [5.0; 1.0]
    @test parentindices(get_Fk(b2obj))[1][3] == 4
    @test b2obj.elements[3].f ≈ [3.0; 2.0] && b2obj.elements[3].h == 0

    # fourth point
    insert_flag = add_feasible!(b2obj, OVector([2.5; 1.5], 0), m2d)
    @test insert_flag == :dominates
    @test length(get_Fk(b2obj)) == 3
    @test get_Fk(b2obj)[3].f == b2obj.elements[5].f ≈ [2.5; 1.5]
    @test parentindices(get_Fk(b2obj))[1][3] == 5
    @test b2obj.elements[4].f ≈ [5.0; 1.0] && b2obj.elements[4].h == 0

    # last point
    insert_flag = add_feasible!(b2obj, OVector([-1.0; 0.0], 0), m2d)
    @test insert_flag == :dominates
    @test length(get_Fk(b2obj)) == 1
    @test get_Fk(b2obj)[1].f == b2obj.elements[6].f ≈ [-1.0; 0.0]
    @test parentindices(get_Fk(b2obj))[1][1] == 6
    @test b2obj.elements[5].f ≈ [2.5; 1.5] && b2obj.elements[5].h == 0

    # Check mesh updating deep copy
    refine_frame_size!(b2obj.meshes[6])
    for i in 1:5
        @test b2obj.meshes[i].b_init == [0; 1]
        @test b2obj.meshes[i].bᵏ == [0; 1]
        @test b2obj.meshes[i].aᵏ == [5; 2]
    end
    @test b2obj.meshes[6].b_init == [0; 1]
    @test b2obj.meshes[6].bᵏ == [0; 1]
    @test b2obj.meshes[6].aᵏ == [2; 1]

    @test length(get_Uk(b2obj)) == length(get_Ik(b2obj)) == 0
end

@testset "Feasible insertion for three objectives" begin
    b3obj = Barrier((5, 3), 10)
    m5d = GranularMesh(5, [1.0; 2.0; 4.0; 5.0; 10.0])

    @test_throws(ErrorException, add_feasible!(b3obj, OVector([1.0; -4.0], 0), m5d))

    insert_flag = add_feasible!(b3obj, OVector([2.0; 3.0 ; 1.0], 0), m5d)
    @test insert_flag == :improves
    @test length(get_Fk(b3obj)) == 1
    @test get_Fk(b3obj)[1].f == b3obj.elements[1].f ≈ [2.0; 3.0 ; 1.0]

    # first tentative
    insert_flag = add_feasible!(b3obj, OVector([4.0; 3.0; 1.0], 0), m5d)
    @test insert_flag == :is_dominated
    @test length(get_Fk(b3obj)) == 1
    @test get_Fk(b3obj)[1].f == b3obj.elements[1].f ≈ [2.0; 3.0 ; 1.0]
    @test b3obj.elements[2].f ≈ [4.0; 3.0; 1.0]

    # second tentative
    insert_flag = add_feasible!(b3obj, OVector([3.0; 2.0; 1.0], 0), m5d)
    @test insert_flag == :extends
    @test length(get_Fk(b3obj)) == 2
    @test get_Fk(b3obj)[2].f == b3obj.elements[3].f ≈ [3.0; 2.0 ; 1.0]
    @test parentindices(get_Fk(b3obj))[1][2] == 3
    @test b3obj.elements[2].f ≈ [4.0; 3.0; 1.0]

    # third point
    insert_flag = add_feasible!(b3obj, OVector([5.0; 1.0 ; 0.0], 0), m5d)
    @test insert_flag == :extends
    @test length(get_Fk(b3obj)) == 3
    @test get_Fk(b3obj)[3].f == b3obj.elements[4].f ≈ [5.0; 1.0 ; 0.0]
    @test parentindices(get_Fk(b3obj))[1][3] == 4
    @test b3obj.elements[3].f ≈ [3.0; 2.0; 1.0]

    # last point
    insert_flag = add_feasible!(b3obj, OVector([-10; 0.0; -4.0], 0), m5d)
    @test insert_flag == :dominates
    @test length(get_Fk(b3obj)) == 1
    @test get_Fk(b3obj)[1].f == b3obj.elements[5].f ≈ [-10; 0.0; -4.0]
    @test parentindices(get_Fk(b3obj))[1][1] == 5
    @test b3obj.elements[4].f ≈ [5.0; 1.0; 0.0]

    # Check mesh deepcopy insertion
    for i in 1:5
        if mod(i, 2) == 1
            refine_frame_size!(b3obj.meshes[i])
        end
    end
    for i in 1:5
        if mod(i, 2) == 1
            @test b3obj.meshes[i].bᵏ == [-1, 0, 0, 0, 0]
            @test b3obj.meshes[i].aᵏ == [5, 1, 2, 2, 5]
        else
            @test b3obj.meshes[i].bᵏ == [0, 0, 0, 0, 1]
            @test b3obj.meshes[i].aᵏ == [1, 2, 5, 5, 1]
        end
    end

    @test length(get_Uk(b3obj)) == length(get_Ik(b3obj)) == 0
end

@testset "Infeasible insertion for two objectives" begin

    m3d = GranularMesh(3, [1.0; 2.0; 10.0])
    b2obj = Barrier((3,2), 10)

    @test_throws(ErrorException, add_infeasible!(b2obj, OVector([1.0; -4.0; 6.0], 3.2), m3d))
    @test_throws(ErrorException, add_infeasible!(b2obj, OVector([1.0; -4.0], 0), m3d))

    b2obj.h_max = 4.5

    # first tentative
    insert_flag = add_infeasible!(b2obj, OVector([2.0; 3.0], 7.6), m3d)
    @test insert_flag == :is_rejected
    @test length(get_Fk(b2obj)) == length(get_Uk(b2obj)) == length(get_Ik(b2obj)) == 0
    @test b2obj.elements[1].f ≈ [2.0; 3.0] && b2obj.elements[1].h ≈ 7.6

    insert_flag = add_infeasible!(b2obj, OVector([2.0; 3.0], 3.5), m3d)
    @test insert_flag == :improves
    @test length(get_Uk(b2obj)) == 1
    @test get_Uk(b2obj)[1].f == get_Ik(b2obj)[1].f == b2obj.elements[2].f ≈ [2.0; 3.0]
    @test get_Uk(b2obj)[1].h == get_Ik(b2obj)[1].h == b2obj.elements[2].h ≈ 3.5
    @test parentindices(get_Uk(b2obj))[1][1] == parentindices(get_Ik(b2obj))[1][1] == 2

    # first tentative
    insert_flag = add_infeasible!(b2obj, OVector([3.0; 3.0], 2.1), m3d)
    @test insert_flag == :no_improvement
    @test length(get_Uk(b2obj)) == 2 && length(get_Ik(b2obj)) == 1
    @test get_Uk(b2obj)[2].f == b2obj.elements[3].f ≈ [3.0; 3.0]
    @test get_Uk(b2obj)[2].h == b2obj.elements[3].h ≈ 2.1
    @test parentindices(get_Uk(b2obj))[1][2] == 3
    @test get_Uk(b2obj)[1].f == get_Ik(b2obj)[1].f == b2obj.elements[2].f ≈ [2.0; 3.0]
    @test get_Uk(b2obj)[1].h == get_Ik(b2obj)[1].h == b2obj.elements[2].h ≈ 3.5
    @test parentindices(get_Uk(b2obj))[1][1] == parentindices(get_Ik(b2obj))[1][1] == 2

    # first bis tentative
    insert_flag = add_infeasible!(b2obj, OVector([3.0; 3.0], 3.7), m3d)
    @test insert_flag == :is_dominated
    @test b2obj.elements[4].f ≈ [3.0; 3.0] && b2obj.elements[4].h ≈ 3.7
    @test length(get_Uk(b2obj)) == 2 && length(get_Ik(b2obj)) == 1
    @test get_Uk(b2obj)[2].f == b2obj.elements[3].f ≈ [3.0; 3.0]
    @test get_Uk(b2obj)[2].h == b2obj.elements[3].h ≈ 2.1
    @test parentindices(get_Uk(b2obj))[1][2] == 3
    @test get_Uk(b2obj)[1].f == get_Ik(b2obj)[1].f == b2obj.elements[2].f ≈ [2.0; 3.0]
    @test get_Uk(b2obj)[1].h == get_Ik(b2obj)[1].h == b2obj.elements[2].h ≈ 3.5
    @test parentindices(get_Uk(b2obj))[1][1] == parentindices(get_Ik(b2obj))[1][1] == 2

    # second tentative
    insert_flag = add_infeasible!(b2obj, OVector([3.0; 2.0], 3.9), m3d)
    @test insert_flag == :extends
    @test length(get_Uk(b2obj)) == 3 && length(get_Ik(b2obj)) == 2
    @test get_Uk(b2obj)[3].f == get_Ik(b2obj)[2].f == b2obj.elements[5].f ≈ [3.0; 2.0]
    @test get_Uk(b2obj)[3].h == get_Ik(b2obj)[2].h == b2obj.elements[5].h ≈ 3.9
    @test parentindices(get_Uk(b2obj))[1][3] == parentindices(get_Ik(b2obj))[1][2] == 5
    @test get_Uk(b2obj)[1].f == get_Ik(b2obj)[1].f == b2obj.elements[2].f ≈ [2.0; 3.0]
    @test get_Uk(b2obj)[1].h == get_Ik(b2obj)[1].h == b2obj.elements[2].h ≈ 3.5
    @test parentindices(get_Uk(b2obj))[1][1] == parentindices(get_Ik(b2obj))[1][1] == 2
    @test get_Uk(b2obj)[2].f == b2obj.elements[3].f ≈ [3.0; 3.0]
    @test get_Uk(b2obj)[2].h == b2obj.elements[3].h ≈ 2.1

    # third point
    insert_flag = add_infeasible!(b2obj, OVector([5.0; 1.0], 1.3), m3d)
    @test insert_flag == :extends
    @test length(get_Uk(b2obj)) == 4 && length(get_Ik(b2obj)) == 3
    @test get_Uk(b2obj)[4].f == get_Ik(b2obj)[3].f == b2obj.elements[6].f ≈ [5.0; 1.0]
    @test get_Uk(b2obj)[4].h == get_Ik(b2obj)[3].h == b2obj.elements[6].h ≈ 1.3
    @test parentindices(get_Uk(b2obj))[1][4] == parentindices(get_Ik(b2obj))[1][3] == 6
    @test get_Uk(b2obj)[3].f == get_Ik(b2obj)[2].f == b2obj.elements[5].f ≈ [3.0; 2.0]
    @test get_Uk(b2obj)[3].h == get_Ik(b2obj)[2].h == b2obj.elements[5].h ≈ 3.9
    @test parentindices(get_Uk(b2obj))[1][3] == parentindices(get_Ik(b2obj))[1][2] == 5

    # fourth point
    insert_flag = add_infeasible!(b2obj, OVector([2.5; 1.5], 2.2), m3d)
    @test insert_flag == :dominates
    @test length(get_Uk(b2obj)) == 4 && length(get_Ik(b2obj)) == 3
    @test get_Uk(b2obj)[4].f == get_Ik(b2obj)[3].f == b2obj.elements[7].f ≈ [2.5; 1.5]
    @test get_Uk(b2obj)[4].h == get_Ik(b2obj)[3].h == b2obj.elements[7].h ≈ 2.2
    @test parentindices(get_Uk(b2obj))[1][4] == parentindices(get_Ik(b2obj))[1][3] == 7
    @test get_Uk(b2obj)[3].f == get_Ik(b2obj)[2].f == b2obj.elements[6].f ≈ [5.0; 1.0]
    @test get_Uk(b2obj)[3].h == get_Ik(b2obj)[2].h == b2obj.elements[6].h ≈ 1.3
    @test parentindices(get_Uk(b2obj))[1][3] == parentindices(get_Ik(b2obj))[1][2] == 6

    # last point
    insert_flag = add_infeasible!(b2obj, OVector([-1.0; 0.0], 0.6), m3d)
    @test insert_flag == :dominates
    @test length(get_Uk(b2obj)) == 1 && length(get_Ik(b2obj)) == 1
    @test get_Uk(b2obj)[1].f == get_Ik(b2obj)[1].f == b2obj.elements[8].f ≈ [-1.0; 0.0]
    @test get_Uk(b2obj)[1].h == get_Ik(b2obj)[1].h == b2obj.elements[8].h ≈ 0.6
    @test parentindices(get_Uk(b2obj))[1][1] == parentindices(get_Ik(b2obj))[1][1] == 8
    @test b2obj.elements[7].f ≈ [2.5; 1.5] && b2obj.elements[7].h ≈ 2.2

    # check deep mesh copy insertion
    refine_frame_size!(m3d)
    for i in 1:7
        @test b2obj.meshes[i].aᵏ != m3d.aᵏ && b2obj.meshes[i].bᵏ != m3d.bᵏ
    end
    for i in 1:7
        if mod(i, 2) == 0
            refine_frame_size!(b2obj.meshes[i])
        end
    end
    for i in 1:7
        if mod(i, 2) == 0
            @test b2obj.meshes[i].aᵏ == m3d.aᵏ && b2obj.meshes[i].bᵏ == m3d.bᵏ
        end
    end

    @test length(get_Fk(b2obj)) == 0
end

@testset "Update barrier" begin

    m2d = GranularMesh(2, [6.0, 1.8])
    b3obj = Barrier((2, 3), 10)

    @test_throws(ErrorException, add_infeasible!(b3obj, OVector([1.0; -4.0], 2.5), m2d))

    insert_flag = add_infeasible!(b3obj, OVector([2.0; 3.0 ; 1.0], 2.5), m2d)
    @test insert_flag == :improves

    # first tentative
    insert_flag = add_infeasible!(b3obj, OVector([4.0; 3.0; 1.0], 4.3), m2d)
    @test insert_flag == :is_dominated

    # second tentative
    insert_flag = add_infeasible!(b3obj, OVector([3.0; 2.0; 1.0], 3.3), m2d)
    @test insert_flag == :extends

    # third point
    insert_flag = add_infeasible!(b3obj, OVector([5.0; 1.0 ; 0.0], 8.7), m2d)
    @test insert_flag == :extends
    @test length(get_Uk(b3obj)) == length(get_Ik(b3obj)) == 3
    @test get_Uk(b3obj)[3].f == get_Ik(b3obj)[3].f == b3obj.elements[4].f ≈ [5.0; 1.0; 0.0]
    @test get_Uk(b3obj)[3].h == get_Ik(b3obj)[3].h == b3obj.elements[4].h ≈ 8.7
    @test parentindices(get_Uk(b3obj))[1][3] == parentindices(get_Ik(b3obj))[1][3] == 4
    @test get_Uk(b3obj)[2].f == get_Ik(b3obj)[2].f == b3obj.elements[3].f ≈ [3.0; 2.0; 1.0]
    @test get_Uk(b3obj)[2].h == get_Ik(b3obj)[2].h == b3obj.elements[3].h ≈ 3.3
    @test parentindices(get_Uk(b3obj))[1][2] == parentindices(get_Ik(b3obj))[1][2] == 3
    @test b3obj.h_max == Inf

    @test_throws(ErrorException, update_barrier!(b3obj, -2.5))

    # first update
    update_barrier!(b3obj, 3.0)
    @test length(get_Uk(b3obj)) == length(get_Ik(b3obj)) == 1
    @test get_Uk(b3obj)[1].f == get_Ik(b3obj)[1].f == b3obj.elements[1].f ≈ [2.0; 3.0; 1.0]
    @test get_Uk(b3obj)[1].h == get_Ik(b3obj)[1].h == b3obj.elements[1].h ≈ 2.5
    @test b3obj.last_index == 4
    @test b3obj.h_max == 3.0

    # another point
    insert_flag = add_infeasible!(b3obj, OVector([6.0; 7.0; 3.0], 0.8), m2d)
    @test insert_flag == :no_improvement
    @test length(get_Uk(b3obj)) == 2 && length(get_Ik(b3obj)) == 1
    @test get_Uk(b3obj)[2].f == b3obj.elements[5].f ≈ [6.0; 7.0; 3.0]
    @test get_Uk(b3obj)[2].h == b3obj.elements[5].h ≈ 0.8
    @test parentindices(get_Uk(b3obj))[1][2] == 5
    @test get_Uk(b3obj)[1].f == get_Ik(b3obj)[1].f == b3obj.elements[1].f ≈ [2.0; 3.0; 1.0]
    @test get_Uk(b3obj)[1].h == get_Ik(b3obj)[1].h == b3obj.elements[1].h ≈ 2.5
    @test parentindices(get_Uk(b3obj))[1][1] == parentindices(get_Ik(b3obj))[1][1] == 1
    @test b3obj.h_max == 3.0

    # and one more again
    insert_flag = add_infeasible!(b3obj, OVector([4.0; 5.0; 2.0], 1.0), m2d)
    @test insert_flag == :no_improvement
    @test length(get_Uk(b3obj)) == 3 && length(get_Ik(b3obj)) == 1
    @test get_Uk(b3obj)[3].f == b3obj.elements[6].f ≈ [4.0; 5.0; 2.0]
    @test get_Uk(b3obj)[3].h == b3obj.elements[6].h ≈ 1.0
    @test parentindices(get_Uk(b3obj))[1][3] == 6
    @test get_Uk(b3obj)[1].f == get_Ik(b3obj)[1].f == b3obj.elements[1].f ≈ [2.0; 3.0; 1.0]
    @test get_Uk(b3obj)[1].h == get_Ik(b3obj)[1].h == b3obj.elements[1].h ≈ 2.5
    @test b3obj.h_max == 3.0

    # second update
    update_barrier!(b3obj, 1.5)
    @test length(get_Uk(b3obj)) == 2
    @test length(get_Ik(b3obj)) == 1
    @test get_Ik(b3obj)[1].f == b3obj.elements[6].f ≈ [4.0; 5.0; 2.0]
    @test get_Ik(b3obj)[1].h == b3obj.elements[6].h ≈ 1.0
    @test parentindices(get_Ik(b3obj))[1][1] == 6
    @test b3obj.h_max == 1.5

    # third update
    update_barrier!(b3obj, 0.9)
    @test length(get_Uk(b3obj)) == 1
    @test length(get_Ik(b3obj)) == 1
    @test get_Uk(b3obj)[1].f == get_Ik(b3obj)[1].f == b3obj.elements[5].f ≈ [6.0; 7.0; 3.0]
    @test get_Ik(b3obj)[1].h == get_Ik(b3obj)[1].h == b3obj.elements[5].h ≈ 0.8
    @test parentindices(get_Uk(b3obj))[1][1] == parentindices(get_Ik(b3obj))[1][1] == 5
    @test b3obj.h_max == 0.9

    # last point
    insert_flag = add_infeasible!(b3obj, OVector([10; 8.0; 4.0], 0.7), m2d)
    @test insert_flag == :no_improvement
    @test length(get_Uk(b3obj)) == 2 && length(get_Ik(b3obj)) == 1
    @test get_Uk(b3obj)[2].f == b3obj.elements[7].f ≈ [10.0; 8.0; 4.0]
    @test get_Uk(b3obj)[2].h == b3obj.elements[7].h ≈ 0.7
    @test parentindices(get_Uk(b3obj))[1][2] == 7
    @test get_Uk(b3obj)[1].f == get_Ik(b3obj)[1].f == b3obj.elements[5].f ≈ [6.0; 7.0; 3.0]
    @test get_Ik(b3obj)[1].h == get_Ik(b3obj)[1].h == b3obj.elements[5].h ≈ 0.8
    @test parentindices(get_Uk(b3obj))[1][1] == parentindices(get_Ik(b3obj))[1][1] == 5

    update_barrier!(b3obj, 0.5)
    @test length(get_Uk(b3obj)) == length(get_Ik(b3obj)) == 0

    update_barrier!(b3obj, 0.0)
    @test length(get_Uk(b3obj)) == length(get_Ik(b3obj)) == 0

    @test length(get_Fk(b3obj)) == 0
end

@testset "Poll center selection for two objectives" begin

    b2obj = Barrier((2, 2), 10)
    m2d = GranularMesh(2, [3.5; 17.5])

    # Build a feasible Pareto front with five elements
    current_frame_centers = frame_centers(b2obj, 0)
    @test current_frame_centers.feasible == 0
    @test current_frame_centers.infeasible == 0

    add_feasible!(b2obj, OVector([1.0; 3.5], 0), m2d)

    current_frame_centers = frame_centers(b2obj, 0)
    @test current_frame_centers.feasible == 1
    @test current_frame_centers.infeasible == 0

    add_feasible!(b2obj, OVector([3.4; 1.5], 0), m2d)

    current_frame_centers = frame_centers(b2obj, 0)
    @test current_frame_centers.feasible == 1
    @test current_frame_centers.infeasible == 0

    add_feasible!(b2obj, OVector([1.5; 3.0], 0), m2d)
    add_feasible!(b2obj, OVector([5.5; 1.0], 0), m2d)
    add_feasible!(b2obj, OVector([6.0; 0.5], 0), m2d)

    current_frame_centers = frame_centers(b2obj, 0)
    @test current_frame_centers.feasible == 2
    @test current_frame_centers.infeasible == 0

    # Modifications of some information
    for i in 1:4
        refine_frame_size!(b2obj.meshes[i])
    end

    current_frame_centers = frame_centers(b2obj, 0)
    @test current_frame_centers.feasible == 5
    @test current_frame_centers.infeasible == 0

    # New modifications
    refine_frame_size!(b2obj.meshes[2])
    refine_frame_size!(b2obj.meshes[5])
    refine_frame_size!(b2obj.meshes[5])

    current_frame_centers = frame_centers(b2obj, 0)
    @test current_frame_centers.feasible == 3
    @test current_frame_centers.infeasible == 0

    # Increase the selection parameter value
    current_frame_centers = frame_centers(b2obj, 4)
    @test current_frame_centers.feasible == 2
    @test current_frame_centers.infeasible == 0

    # Reset the information
    for i in 1:5
        b2obj.meshes[i] = deepcopy(m2d)
    end

    # Check infeasibility selection
    add_infeasible!(b2obj, OVector([4.0;1.0], 15.4), m2d)

    current_frame_centers = frame_centers(b2obj, 0)
    @test current_frame_centers.feasible == 2
    @test current_frame_centers.infeasible == 6

    add_infeasible!(b2obj, OVector([3.0;4.0], 16.8), m2d)

    current_frame_centers = frame_centers(b2obj, 0)
    @test current_frame_centers.feasible == 2
    @test current_frame_centers.infeasible == 6

    add_infeasible!(b2obj, OVector([3.2;1.3], 14.4), m2d)
    current_frame_centers = frame_centers(b2obj, 0; use_doM_selection=false)
    @test current_frame_centers.feasible == 2
    @test current_frame_centers.infeasible == 8

    current_frame_centers = frame_centers(b2obj, 0)
    @test current_frame_centers.feasible == 2
    @test current_frame_centers.infeasible == 6

    # Remove articially feasible points and check infeasibility selection
    b2obj.within_Fk[1:b2obj.last_index] .= false
    current_frame_centers = frame_centers(b2obj, 0)
    @test current_frame_centers.feasible == 0
    @test current_frame_centers.infeasible == 7

    refine_frame_size!(b2obj.meshes[7])
    current_frame_centers = frame_centers(b2obj, 0)
    @test current_frame_centers.feasible == 0
    @test current_frame_centers.infeasible == 6

    refine_frame_size!(b2obj.meshes[6])
    current_frame_centers = frame_centers(b2obj, 0)
    @test current_frame_centers.feasible == 0
    @test current_frame_centers.infeasible == 8
    @test length(get_Ik(b2obj)) == length(get_Uk(b2obj)) == 3
end

@testset "Saving barrier values" begin
    b3obj = Barrier((2,3), 10)
    m2d = GranularMesh(2, [3.5; 17.5])

    @test add_feasible!(b3obj, OVector([2.0; 3.0; 1.0], 0), m2d) == :improves
    @test add_feasible!(b3obj, OVector([3.0; 2.0; 1.0], 0), m2d) == :extends
    @test add_feasible!(b3obj, OVector([5.0; 1.0; 0.0], 0), m2d) == :extends

    save_pf_values(b3obj ,"simple-barrier-test.txt")
    counter = 0
    open("simple-barrier-test.txt", "r") do io
        for ln in eachline(io)
            counter += 1
        end
    end
    @test counter == 3
    rm("simple-barrier-test.txt")

end
