using DMultiMadsPB

# Define functions and constraints
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

# Create a blackbox problem
constrained_Kursawe = BBProblem(x -> [Kursawe_obj(x); constraints(x)],
                                3, 4, # number of inputs, number of outputs
                                [OBJ, OBJ, CSTR, CSTR], # type of outputs
                                lvar=-5 * ones(3), # lower bounds
                                uvar = 5 * ones(3), # upper bounds
                                name="Kursawe_constrained")

# Create a Mads model
model = MadsModel(constrained_Kursawe)

# You can fix options (which are already configured by default)
model.options.neval_bb_max = 2500 # Maximum budget allowed for the solver
# model.options.display = false # Display the optimization process (true by default)
# model.options.is_ordered = false # Order candidates
# model.options.is_opportunistic = false

# And some other (take a look at the module)

# By default, the DMulti-MADS-PB strategy is used. To change,
# model.options.use_penalty_approach = true # To active penalty approach
# model.params.h_init = 0 # To activate two phase approach

# Choose start points
start_points = [constrained_Kursawe.meta.lvar[:,] +  (j - 1) * (constrained_Kursawe.meta.uvar[:,] - constrained_Kursawe.meta.lvar[:,]) / (constrained_Kursawe.meta.ninputs - 1)  for j in 1:constrained_Kursawe.meta.ninputs]

# Solve problem
stop = solve!(model, start_points)

# You can get more information after the resolution

# for example the Cache
println(model.cache)

# Or the best feasible solutions
println(model.cache.Vᵏ[parentindices(get_Fk(model.barrier))[1]])
