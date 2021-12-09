@testset "Mads algorithm : Rosenbrock function" begin

    Rosenbrock = BBProblem(x -> [(x[1] - 1)^2 + 100 * (x[2] - x[1]^2)^2],
                           2, 1, [OBJ], lvar=[-5.0; -5.0], uvar=[5.0; 5.0],
                           name="Rosenbrock")

    @test Rosenbrock.meta.name == "Rosenbrock"
    @test Rosenbrock.meta.typeoutputs == [OBJ]
    @test Rosenbrock.meta.ninputs == 2
    @test Rosenbrock.meta.noutputs == 1
    @test Rosenbrock.meta.lvar ≈ [-5.0; -5.0]
    @test Rosenbrock.meta.uvar ≈ [5.0; 5.0]
    @test eval_x(Rosenbrock, [1.0; 1.0]) ≈ [0]

    model = MadsModel(Rosenbrock)
    model.options.neval_bb_max = 1000

    stop = solve!(model, [[-1.0; 4.0]])
    @test stop == MIN_MESH_REACHED
    @test model.cache.last_index == 630
    @test length(get_Fk(model.barrier)) == 1
    @test abs(get_Fk(model.barrier)[1].f[1]) <= 1e-8
    sol = model.cache.Vᵏ[parentindices(get_Fk(model.barrier))[1][1]]
    @test abs(sol.inputs[1] - sol.inputs[2]) <= 1e-4
end

@testset "Mads algorithm: the mustache example (proposed for PB)" begin
    # objectives
    function f(x)
        return sqrt((x[1]-20)^2 + (x[2]-1)^2)
    end

    # Constraints
    function c(x)
        return [sin(x[1]) - 1/10 - x[2], x[2] - sin(x[1])]
    end

    Mustache = BBProblem(x -> [f(x);c(x)],
                         2, 3, [OBJ; CSTR; CSTR],
                         lvar=[-2.0;-2.0], uvar=[22.0;2.0],
                         name = "Mustache")

    @test Mustache.meta.name == "Mustache"
    @test Mustache.meta.typeoutputs == [OBJ; CSTR; CSTR]
    @test Mustache.meta.ninputs == 2
    @test Mustache.meta.noutputs == 3
    @test Mustache.meta.lvar ≈ [-2.0;-2.0]
    @test Mustache.meta.uvar ≈ [22.0;2.0]

    for h_init in [Inf, 0]
        model = MadsModel(Mustache)
        model.options.neval_bb_max = 2000
        model.options.is_opportunistic = false
        model.params.h_init = h_init

        stop = solve!(model, [[0.0;0.0]])
        @test isapprox(norm(get_Fk(model.barrier)[1].f, Inf), 0.08, atol=1e-2)
    end
end

@testset "Mads algorithm : a simple example 2 objectives" begin

    SP1 = BBProblem(x -> [(x[1] - 1)^2 + (x[1] - x[2])^2;
                          (x[1] - x[2])^2 + (x[2] - 3)^2],
                    2, 2, [OBJ; OBJ], lvar=[-1.0; -1.0], uvar=[5.0; 5.0],
                    name="SP1")

    @test SP1.meta.name == "SP1"
    @test SP1.meta.typeoutputs == [OBJ; OBJ]
    @test SP1.meta.ninputs == 2
    @test SP1.meta.noutputs == 2
    @test SP1.meta.lvar ≈ [-1.0; -1.0]
    @test SP1.meta.uvar ≈ [5.0; 5.0]
    @test eval_x(SP1, [1.5; 1.5]) ≈ [0.25; 2.25]

    println("Running two times DMultiMads...")
    for penalty_approach in [true, false]
        model = MadsModel(SP1)
        model.options.neval_bb_max = 2000
        model.options.display = false
        model.options.use_penalty_approach = penalty_approach
        model.params.w_min = 0
        stop = solve!(model, [[1.5; 1.5]])
        @test stop == MAX_BB_REACHED
        @test length(get_Fk(model.barrier)) == 259
    end
end

@testset "Mads algorithm : a simple example 3 objectives" begin

    Viennet = BBProblem(x -> [0.5 * (x[1]^2 + x[2]^2) + sin(x[1]^2 + x[2]^2);
                              (3 * x[1] - 2 * x[2] + 4)^2 / 8 + (x[1] - x[2] + 1)^2 / 27 + 15;
                              1 / (x[1]^2 + x[2]^2 + 1) - 1.1 * exp(-(x[1]^2 + x[2]^2))],
                        2, 3, [OBJ; OBJ; OBJ], lvar=[-3.0; -3.0], uvar=[3.0; 3.0],
                        name="Viennet")

    @test Viennet.meta.name == "Viennet"
    @test Viennet.meta.typeoutputs == [OBJ; OBJ; OBJ]
    @test Viennet.meta.ninputs == 2
    @test Viennet.meta.noutputs == 3
    @test Viennet.meta.lvar ≈ [-3.0; -3.0]
    @test Viennet.meta.uvar ≈ [3.0; 3.0]
    @test eval_x(Viennet, [1.0; 1.0]) ≈ [1 + sin(2); 25 / 8 + 1 / 27 + 15; 1/3 - 1.1  * exp(-2)]

    println("Running model two times with different variants")
    display = true
    for penalty_approach in [true, false]
        model = MadsModel(Viennet)
        model.options.neval_bb_max = 2000
        model.options.use_penalty_approach = penalty_approach
        model.options.is_opportunistic = false
        model.options.is_ordered = false
        model.options.display = display
        display = !display
        stop = solve!(model, [[1.5; 1.5]])
        @test stop == MAX_BB_REACHED
        @test length(model.cache.Vᵏ[1:model.cache.last_index]) == 2000
        @test length(get_Fk(model.barrier)) == 539
    end
end

@testset "Mads algorithm: FES1 with constraints" begin

    function constraints(x)
        n = length(x)
        return (3 .- 2 * x[2:n-1]) .* x[2:n-1] - x[1:n-2] - 2 * x[3:n] .+ 1
    end

    function FES1_obj(x)
        n = 10;

        # objective funct(ion
        f1 = sum(abs.(x - exp.(((1:n)/n).^2) /3).^0.5);
        f2 = sum((x - 0.5 * cos.(10 * pi * (1:n)/ n) .- 0.5).^2);

        return [f1;f2]
    end

    constrained_FES1 = BBProblem(x -> [FES1_obj(x); constraints(x)],
                                 10, 10, [[OBJ; OBJ]; repeat([CSTR], 8)],
                                 lvar=zeros(10), uvar = ones(10),
                                 name="FES1_constrained")
    @test constrained_FES1.meta.name == "FES1_constrained"
    @test constrained_FES1.meta.typeoutputs == [OBJ; OBJ; CSTR; CSTR; CSTR; CSTR; CSTR; CSTR ; CSTR; CSTR]
    @test constrained_FES1.meta.ninputs == 10
    @test constrained_FES1.meta.noutputs == 10
    @test constrained_FES1.meta.lvar ≈ zeros(10)
    @test constrained_FES1.meta.uvar ≈ ones(10)

    nb_solutions_found = Int[]
    for h_init in [Inf, 0]
        model = MadsModel(constrained_FES1)
        model.options.neval_bb_max = 5000
        model.params.h_init = h_init
        model.options.display = false
        stop = solve!(model, [(constrained_FES1.meta.lvar + constrained_FES1.meta.uvar) / 2])
        push!(nb_solutions_found, length(get_Fk(model.barrier)))
        @test stop == MAX_BB_REACHED
    end
    @test nb_solutions_found[1] > nb_solutions_found[2]
end

@testset "Mads algorithm: Kursawe with constraints 4" begin

    function constraints(x)
        n = length(x)
        return x[1:n-1].^2 + x[2:n].^2 + x[1:n-1] .* x[2:n] .- 1
    end

    function Kursawe_obj(x)
        #params
        n = 3;

        # functions
        f1 = sum(-10 * exp.(-0.2 * sqrt.(x[1:n-1].^2 + x[2:n].^2)));
        f2 = sum(abs.(x[1:n]).^0.8 + 5 * sin.(x[1:n]).^3);

        return [f1;f2]

    end

    constrained_Kursawe = BBProblem(x -> [Kursawe_obj(x); constraints(x)],
                                    3, 4, [OBJ, OBJ, CSTR, CSTR],
                                    lvar=-5 * ones(3), uvar = 5 * ones(3),
                                    name="Kursawe_constrained")
    @test constrained_Kursawe.meta.name == "Kursawe_constrained"
    @test constrained_Kursawe.meta.typeoutputs == [OBJ, OBJ, CSTR, CSTR]
    @test constrained_Kursawe.meta.ninputs == 3
    @test constrained_Kursawe.meta.noutputs == 4
    @test constrained_Kursawe.meta.lvar ≈ -5 * ones(3)
    @test constrained_Kursawe.meta.uvar ≈ 5 * ones(3)
    model = MadsModel(constrained_Kursawe)
    model.options.neval_bb_max = 2500
    model.options.display = true
    stop = solve!(model, [(constrained_Kursawe.meta.lvar + constrained_Kursawe.meta.uvar) / 2])
    @test stop == MAX_BB_REACHED
    @test length(get_Fk(model.barrier)) ==  182
    @test length(get_Ik(model.barrier)) ==  1
end

@testset "Mads algorithm: a simple example with constraints" begin

    function constraints(x)
        n = length(x)
        return (3 .- 2 * x[2:n-1]) .* x[2:n-1] - x[1:n-2] - 2 * x[3:n] .+ 1
    end

    function CL1_obj(x)
        F = 10
        E = 2*10^5
        L = 200
        sigma = 10
        return [L * (2 * x[1] + sqrt(2) * x[2] + sqrt(x[3]) + x[4]);
                F * L / E * (2 / x[1] + (2 * sqrt(2)) / x[2] - (2 * sqrt(2)) / x[3] + 2 / x[4])]
    end

    F_p = 10.0
    E_p = 2.0 * 10^5
    L_p = 200.0
    σ_p = 10

    constrained_CL1 = BBProblem(x -> [CL1_obj(x); constraints(x)],
                                4, 4, [OBJ; OBJ; CSTR; CSTR],
                                lvar=[F_p / σ_p ; sqrt(2) * F_p / σ_p;
                                      sqrt(2) * F_p / σ_p ; F_p / σ_p],
    uvar = (3 * F_p / σ_p) * ones(4),
    name="CL1_constrained")

    @test constrained_CL1.meta.name == "CL1_constrained"
    @test constrained_CL1.meta.typeoutputs == [OBJ; OBJ; CSTR; CSTR]
    @test constrained_CL1.meta.ninputs == 4
    @test constrained_CL1.meta.noutputs == 4
    @test constrained_CL1.meta.lvar ≈ [F_p / σ_p ; sqrt(2) * F_p / σ_p;
                                       sqrt(2) * F_p / σ_p ; F_p / σ_p]
    @test constrained_CL1.meta.uvar ≈ (3 * F_p / σ_p) * ones(4)

    model = MadsModel(constrained_CL1)
    model.options.neval_bb_max = 2500

    stop = solve!(model, [(constrained_CL1.meta.lvar + constrained_CL1.meta.uvar) / 2])
    @test stop == MAX_BB_REACHED
    @test length(get_Fk(model.barrier)) == 517
end
