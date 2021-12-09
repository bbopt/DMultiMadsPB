export MadsModel, solve!, add_step!, clear_steps!
export UNDEFINED, NO_INIT_CANDIDATES, STOP_IF_FEASIBLE, MIN_MESH_REACHED, MAX_BB_REACHED, MAX_BB_OUTBOUND_REACHED
export NO_SUCCESS, PARTIAL_SUCCESS, FULL_SUCCESS

# Reasons for the algorithm to stop
@enum StopReason begin
    UNDEFINED
    NO_INIT_CANDIDATES
    STOP_IF_FEASIBLE
    MIN_MESH_REACHED
    MAX_BB_REACHED
    MAX_BB_OUTBOUND_REACHED
end

# Success type at the end of the iteration
@enum SuccessType NO_SUCCESS PARTIAL_SUCCESS FULL_SUCCESS

######################################################
# Private main constants
######################################################
const available_poll_symbols = [:poll_2nd, :poll_np1d, :poll_2d, :poll_1d]
const available_search_symbols = [:search_speculative]
const symbols_to_stepsets = (poll_2nd = PollSet2N,
                             poll_np1d = PollSetNp1,
                             poll_2d = PollSet2,
                             poll_1d = PollSet1,
                             search_speculative = SpeculativeSearchSet)

######################################################
# Configuration options of a Mads model
######################################################

# Model statistics
mutable struct MadsStatistics
    neval_bb :: Int # number of calls to the blackbox
    niterations:: Int # number of iterations
    ncache_hits :: Int # number of cache hits
    noutbound_hits :: Int # number of points evaluated outbounds

    neval_bb_feasible :: Int # number of points evaluated which are feasible
    neval_bb_infeasible :: Int # number of points evaluated which are infeasible

    nfull_successes::Int # number of successful iterations
    npartial_successes::Int # number of partial iterations
    nno_successes::Int # number of failed iterations
    nopportunistic_triggers::Int # number of opportunistic iterations

    function MadsStatistics()
        return new(0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
    end
end

# Model internal parameters
mutable struct MadsParameters

    h_init :: Float64
    ρ_trigger :: Float64
    max_size :: Int
    w_min:: Int
    seed :: Int
    rng :: AbstractRNG

    function MadsParameters(;
                            h_init::Float64=Inf,
                            ρ_trigger::Float64=0.1,
                            max_size::Int= 30000,
                            w_min::Int=1,
                            seed::Int=1234)
        return new(h_init, ρ_trigger, max_size, w_min, seed, MersenneTwister(seed))
    end

end

function check_mads_parameters!(params::MadsParameters)
    if params.ρ_trigger <= 0
        error("ρ_trigger is a positive parameter")
    end
    if params.w_min < 0
        error("w_min is a positive parameter")
    end
    if params.h_init < 0
        error("h_init is a positive parameter")
    end
    if params.max_size <= 0
        error("max_size is a strictly positive parameter")
    end

    # reset rng in case someone changes the seed
    params.rng = MersenneTwister(params.seed)
end

# Mads iteration attributes
mutable struct MadsIterationAttributes

    # Steps that must be executed around each iteration center
    first_center_steps :: Vector{Symbol}
    second_center_steps :: Vector{Symbol}

    function MadsIterationAttributes()
        return new([:search_speculative, :poll_np1d],
                   #  [:search_speculative, :poll_1d])
                   [:search_speculative, :poll_2d])
    end

end

function check_mads_iteration_attributes(iter_attributes::MadsIterationAttributes)
    if !any(step -> step in available_poll_symbols, iter_attributes.first_center_steps)
        error("A poll step must be performed around the first center")
    end
    if !all(step -> step in union(available_poll_symbols, available_search_symbols),
            iter_attributes.first_center_steps)
        error("Unknown step symbol used in first_center_steps instructions")
    end
    if !all(step -> step in union(available_poll_symbols, available_search_symbols),
            iter_attributes.second_center_steps)
        error("Unknown step symbol used in second_center_steps instructions")
    end
end

mutable struct MadsOptions

    min_tol::Float64 # minimum tolerance
    neval_bb_max::Int # maximum blackbox calls
    noutbound_hits_max::Int # maximum number of evaluations outside bounds

    is_ordered::Bool # order points during the iteration
    is_opportunistic::Bool # opportunistic: stops the iteration as soon as a new point is found

    use_dms_success::Bool # use DMS success condition, similar to DMS, for Fk
    use_Nomad_partial_success::Bool # Nomad does not check after the iteration if there are some points in
                                    # the cache with a h value inferior to the infeasible frame center.
                                    # Practically, it is more efficient, even if it does not correspond to
                                    # the article description (but the convergence remains the same).
    use_doM_trigger::Bool # use dominance move to order frame centers; if set to false, order frame centers
                          # based on Pareto dominance

    use_penalty_approach::Bool # use penalty approach to deal with inequality constraints, similar to DFMO

    display::Bool

    function MadsOptions(;
                         min_tol::Float64=10^(-9),
                         neval_bb_max::Int=20000,
                         noutbound_hits_max::Int=60000,
                         is_ordered::Bool=true,
                         is_opportunistic::Bool=true,
                         use_dms_success::Bool=false,
                         use_Nomad_partial_success::Bool=true,
                         use_doM_trigger::Bool=true,
                         use_penalty_approach::Bool=false,
                         display::Bool=true)
        return new(min_tol, neval_bb_max, noutbound_hits_max, is_ordered,
                   is_opportunistic, use_dms_success, use_Nomad_partial_success,
                   use_doM_trigger, use_penalty_approach,
                   display)
    end

end

function check_mads_options(options::MadsOptions)
    if options.neval_bb_max <= 0
        error("The number of blackbox evaluations cannot be negative !")
    end
    if options.noutbound_hits_max <= 0
        error("The number of blackbox tentative evaluations outside the bounds cannot be negative !")
    end
end


# Keep information relatively to one iteration
mutable struct MadsState

    Fk_frame_center::Int
    Uk_frame_center::Int

    ordered_frame_centers::Tuple{Int, Int}

    h_max::Union{Float64, Nothing}

    is_phase_one::Bool

    last_success::SuccessType
    stop_reason::StopReason

    function MadsState()
        return new(0,0,(0,0),
                   nothing, false,
                   NO_SUCCESS,UNDEFINED)
    end
end

######################################################
# MadsModel
######################################################
mutable struct MadsModel

    bbproblem::BBProblem

    # Parameters and options
    params::MadsParameters
    options::MadsOptions
    attributes::MadsIterationAttributes

    state::MadsState

    # Storing
    barrier::Barrier
    cache::Cache

    # Stats
    stats::MadsStatistics

    function MadsModel(bb::BBProblem)
        # By convention, a MadsModel must solve a problem with upper and lower bounds
        if any(max.(bb.meta.uvar) == Inf) || any(min.(bb.meta.lvar) == -Inf)
            error("Problem does not possess finite lower bounds or upper bounds !")
        end

        n = bb.meta.ninputs
        m = count(==(OBJ), bb.meta.typeoutputs)

        return new(bb,
                   MadsParameters(),
                   MadsOptions(),
                   MadsIterationAttributes(),
                   MadsState(),
                   Barrier((n,m), MadsParameters().max_size; h_max=MadsParameters().h_init),
                   Cache((n, bb.meta.noutputs), MadsParameters().max_size),
                   MadsStatistics())
    end
end

###### Stopping criteria functions ######
function reach_nevals_bb_max(model::MadsModel)
   return (model.stats.neval_bb >= model.options.neval_bb_max)
end

function reach_noutbound_hits_max(model::MadsModel)
    return (model.stats.noutbound_hits >= model.options.noutbound_hits_max)
end


###### Main iteration functions (no exportable) ######
function eval!(model::MadsModel, x::Vector{Float64}, m::GranularMesh, parent_ind::Int)
    if isincache(model.cache, x)
        model.stats.ncache_hits += 1
        return nothing
    end

    if (any(x .> model.bbproblem.meta.uvar) || any(x .< model.bbproblem.meta.lvar))
        println("Warning : candidate out of bounds : ", x)
        model.stats.noutbound_hits += 1
        return nothing
    end

    bb_outputs = eval_x(model.bbproblem, x)
    model.stats.neval_bb += 1

    add_cache!(model.cache, EvalPoint(x, bb_outputs))

    t_outputs = model.bbproblem.meta.typeoutputs
    nb_obj = count(==(OBJ), t_outputs)
    index_obj = findall(==(OBJ), t_outputs)
    index_cstrs = findall(==(CSTR), t_outputs)

    v = begin
        # Follow DFMO implementation of penalty function
        if model.options.use_penalty_approach
            penalty_val = 0
            for ind in index_cstrs
                if max(0, bb_outputs[ind]) < 1
                    penalty_val += 1000 * max(0, bb_outputs[ind])
                else
                    penalty_val += 10 * max(0, bb_outputs[ind])
                end
            end
            OVector(bb_outputs[index_obj] .+ penalty_val, 0)
        else # Filter approach
            OVector(bb_outputs[index_obj],
                    sum(max.(0, bb_outputs[index_cstrs]).^2))
        end
    end

    # The evaluation can return sometimes NaN.
    # In this case, return declare the vector as an Inf.
    if any(isnan, v.f) || isnan(v.h)
        v = OVector(Inf * ones(length(index_obj)), Inf)
    end

    insertion_flag = begin
        if isFeasible(v)
            model.stats.neval_bb_feasible += 1
            add_feasible!(model.barrier, v, m)
        else
            model.stats.neval_bb_infeasible += 1
            add_infeasible!(model.barrier, v, m)
        end
    end

    model.barrier.parent_indexes[model.barrier.last_index] = parent_ind
    return insertion_flag
end

### Initialization phase ###
function init_phase!(model::MadsModel, x0s::Vector{Vector{Float64}})

    # Preliminary checks
    check_mads_parameters!(model.params)
    check_mads_options(model.options)
    check_mads_iteration_attributes(model.attributes)

    # Reset storing structures if change in parameters
    model.params.max_size = max(model.params.max_size, model.options.neval_bb_max)
    model.barrier = Barrier((model.bbproblem.meta.ninputs, count(==(OBJ), model.bbproblem.meta.typeoutputs)),
                            model.params.max_size;
                            h_max = model.params.h_init)
    model.cache = Cache((model.bbproblem.meta.ninputs, model.bbproblem.meta.noutputs),
                        model.params.max_size)

    model.state.stop_reason = NO_INIT_CANDIDATES

    # Evaluate candidates
    for x0 in x0s
        mesh = GranularMesh(model.bbproblem.meta.ninputs,
                            0.1 * (model.bbproblem.meta.uvar - model.bbproblem.meta.lvar))
        result = eval!(model, x0, mesh, 0)

        if reach_nevals_bb_max(model)
            model.state.stop_reason = MAX_BB_REACHED
            break
        end

        if reach_noutbound_hits_max(model)
            model.state.stop_reason = MAX_BB_OUTBOUND_REACHED
            break
        end
    end

    # No start candidates: can occur when x0s is empty
    if model.cache.last_index == 0
        return
    end

    # No init candidates: all evaluated points are above the h_max threshold.
    # Trigger phase one.
    if length(get_Fk(model.barrier)) == 0 && length(get_Uk(model.barrier)) == 0
        model.state.stop_reason = UNDEFINED
        model.state.is_phase_one = true
        return
    end

    # Other cases
    if model.state.stop_reason ∉ [MAX_BB_REACHED, MAX_BB_OUTBOUND_REACHED]
        model.state.stop_reason = UNDEFINED
        return
    end
end

### Search and poll phase functions ###
function search_and_poll!(model::MadsModel)

    # Set mesh
    Mk = begin
        if iszero(model.state.Fk_frame_center)
            model.barrier.meshes[model.state.ordered_frame_centers[1]]
        else
            model.barrier.meshes[model.state.Fk_frame_center]
        end
    end

    if check_mesh_for_stopping(Mk)
        model.state.stop_reason = MIN_MESH_REACHED
        return
    end

    # Generate candidates
    candidates = Array{Float64}(undef, model.bbproblem.meta.ninputs, 0)
    parent_index_candidates = Int[]
    generated_during_search_step = Bool[]
    for frame_center_id in 1:2
        steps = getfield(model.attributes, frame_center_id)
        frame_center = model.state.ordered_frame_centers[frame_center_id]
        for step in steps
            if !iszero(frame_center)
                x_frame = model.cache.Vᵏ[frame_center].inputs
                parent_index = model.barrier.parent_indexes[frame_center]
                last_success_direction = parent_index == 0 ? nothing : x_frame - model.cache.Vᵏ[parent_index].inputs

                is_search_step = false

                generated_candidates_during_step = begin
                    stepset = symbols_to_stepsets[step](Mk, x_frame, last_success_direction=last_success_direction)
                    is_search_step = typeof(stepset) <: AbstractSearchSet
                    if model.options.is_ordered
                        generate_ordered_candidates(stepset, model.params.rng)
                    else
                        generate_candidates(stepset, model.params.rng)
                    end
                end
                # Snap candidates to bounds if possible
                generated_candidates_during_step = begin
                    if isempty(generated_candidates_during_step)
                        generated_candidates_during_step
                    else
                        mapslices(x_candidate -> project_on_mesh_and_snap_to_bounds(Mk, x_candidate, x_frame,
                                                                                    model.bbproblem.meta.lvar,
                                                                                    model.bbproblem.meta.uvar),
                                  generated_candidates_during_step, dims=1)
                    end
                end
                # Add generated points to the list of candidates
                if !isempty(generated_candidates_during_step)
                    candidates = hcat(candidates, generated_candidates_during_step)
                    for i in 1:size(generated_candidates_during_step, 2)
                        push!(parent_index_candidates, frame_center)
                        push!(generated_during_search_step, is_search_step)
                    end
                end
            end
        end
    end

    # Evaluate candidates
    for id in 1:size(candidates, 2)

        insertion_flag = eval!(model, candidates[:, id], Mk,
                               parent_index_candidates[id])
        success_flag = isnothing(insertion_flag) ? NO_SUCCESS : compute_success(model, insertion_flag)

        # Update mesh; slightly different from the article for better performance.
        if !isnothing(insertion_flag) && insertion_flag in [:dominates, :extends]
            x_parent = model.cache.Vᵏ[parent_index_candidates[id]].inputs
            enlarge_frame_size!(model.barrier.meshes[model.barrier.last_index],
                                model.cache.Vᵏ[model.cache.last_index].inputs - x_parent)
        end

        # Update success flag
        if success_flag > model.state.last_success
            model.state.last_success = success_flag
        end

        # Detect potential stopping reasons
        if reach_nevals_bb_max(model)
            model.state.stop_reason = MAX_BB_REACHED
            return
        end
        if reach_noutbound_hits_max(model)
            model.state.stop_reason = MAX_BB_OUTBOUND_REACHED
            return
        end

        # Detect feasibility in case we are in phase one.
        if model.state.is_phase_one && !isnothing(insertion_flag)
            # TODO check if could not be done according to h_start
            if model.barrier.elements[model.barrier.last_index].h == 0
                model.state.stop_reason = STOP_IF_FEASIBLE
                return
            end
        end

        # Trigger opportunistic strategy only when an iteration is considered
        # as a full success
        if (success_flag > PARTIAL_SUCCESS) && model.options.is_opportunistic
            model.stats.nopportunistic_triggers += 1
            break
        end
    end
end

function compute_success(model::MadsModel, insertion_flag)

    v = model.barrier.elements[model.barrier.last_index]

    # Phase one
    if model.state.is_phase_one
        if v.h < model.barrier.elements[model.state.ordered_frame_centers[1]].h
            return FULL_SUCCESS
        else
            return NO_SUCCESS
        end
    else
        success_flag = NO_SUCCESS
        # Feasible case: full success as soon as a new non dominated point which dominates
        # the current feasible frame center is generated.
        if isFeasible(v)
            # The set of feasible points can be empty before the insertion of the new point.
            if !iszero(model.state.Fk_frame_center)
                if model.options.use_dms_success && insertion_flag in [:dominates, :extends, :improves]
                    success_flag = FULL_SUCCESS
                end
                if v <= model.barrier.elements[model.state.Fk_frame_center]
                    success_flag = FULL_SUCCESS
                end
            # If Fk is empty, consider it as a full success, similar to the Nomad software.
            else
                success_flag = FULL_SUCCESS
            end
        else
            # The first iterations are considered as a partial success, when the progressive barrier
            # approach is chosen.
            if isnothing(model.state.h_max)
                if !iszero(model.barrier.h_max)
                    success_flag = PARTIAL_SUCCESS
                end
            else
                if !iszero(model.state.Uk_frame_center)
                    # Partial success if h(x) below h(x_inf)
                    if v.h < model.barrier.elements[model.state.Uk_frame_center].h
                        success_flag = PARTIAL_SUCCESS
                    end
                end
                # Success if change in Iᵏ with h(x) <= h_max.
                if v <= model.barrier.elements[model.state.Uk_frame_center]
                    success_flag = FULL_SUCCESS
                end
                #  if v.h <= model.barrier.elements[model.state.Uk_frame_center].h && insertion_flag in [:dominates, :extends, :improves]
                #      success_flag = FULL_SUCCESS
                #  end
            end
        end
        return success_flag
    end
end

### Update phase functions ###
function set_frame_centers_and_hvalues!(model::MadsModel)
    # First case: phase one has been triggered.
    # There is only one primary frame center, the one with minimum h value.
    if model.state.is_phase_one
        model.state.ordered_frame_centers = (argmin(map(elt -> elt.h, model.barrier.elements[1:model.barrier.last_index])), 0)
    else
        Fk_index, Uk_index = frame_centers(model.barrier, model.params.w_min, use_doM_selection=model.options.use_doM_trigger)

        # Set primary and secondary frame centers according to trigger conditions.
        if iszero(Fk_index)
            model.state.ordered_frame_centers = (Uk_index, 0)
        else
            if iszero(Uk_index)
                model.state.ordered_frame_centers = (Fk_index, 0)
            else
                if model.options.use_doM_trigger
                    # Using extent is slightly more efficient
                    doM = minimum(elt -> sum(elt.f - min.(model.barrier.elements[Uk_index].f, elt.f)), get_Fk(model.barrier))
                    #  if doM >= model.params.ρ_trigger * model.bbproblem.meta.noutputs
                    if doM >= model.params.ρ_trigger * extent(model.barrier)
                        model.state.ordered_frame_centers = (Uk_index, Fk_index)
                    else
                        model.state.ordered_frame_centers = (Fk_index, Uk_index)
                    end
                else # Classic alternative based on dominance but slightly less efficient.
                    if all(model.barrier.elements[Fk_index].f .- model.params.ρ_trigger .>= model.barrier.elements[Uk_index].f)
                        model.state.ordered_frame_centers = (Uk_index, Fk_index)
                    else
                        model.state.ordered_frame_centers = (Fk_index, Uk_index)
                    end
                end
            end
        end

        model.state.Fk_frame_center = Fk_index
        model.state.Uk_frame_center = Uk_index

        # set h_max
        if !isempty(get_Ik(model.barrier))
            model.state.h_max = maximum(map(elt -> elt.h, get_Ik(model.barrier)))
        end
    end
end

# Called at the beginning of each iteration
function update!(model::MadsModel)

    # Phase one is over.
    if model.state.stop_reason == STOP_IF_FEASIBLE
        model.state.is_phase_one = false
        model.state.stop_reason = UNDEFINED
    end

    # No need to update in this case
    if model.state.stop_reason ∈ [MIN_MESH_REACHED, MAX_BB_REACHED, MAX_BB_OUTBOUND_REACHED]
        return
    end

    # There can remain some points in the set Uk which have better h-value,
    # To check if the flag is activated
    if !model.options.use_Nomad_partial_success
        if model.state.last_success == NO_SUCCESS && !iszero(model.state.Uk_frame_center)
            tmp_hx_inf_min = minimum(elt -> elt.h, get_Uk(model.barrier))
            if tmp_hx_inf_min < model.barrier.elements[model.state.Uk_frame_center].h
                model.state.last_success = PARTIAL_SUCCESS
            end
        end
    end

    # Update barrier and mesh
    if model.state.last_success == NO_SUCCESS

        # Set to null the last success directions of the current incumbents.
        for frame_center in 1:2
            barrier_index = model.state.ordered_frame_centers[frame_center]
            if !iszero(barrier_index)
                model.barrier.parent_indexes[barrier_index] = 0
            end
        end

        # Update the mesh
        if !iszero(model.state.Fk_frame_center)
            refine_frame_size!(model.barrier.meshes[model.state.Fk_frame_center])
        else
            refine_frame_size!(model.barrier.meshes[model.state.ordered_frame_centers[1]])
        end

        # Update the barrier threshold
        if !isnothing(model.state.h_max)
            below_hmax_elements = filter(elt -> elt.h < model.state.h_max,
                                         model.barrier.elements[model.barrier.within_Uk .== true])
            if !isempty(below_hmax_elements)
                h_max_tmp = maximum(map(elt-> elt.h, below_hmax_elements))
                if h_max_tmp > model.barrier.elements[model.state.Uk_frame_center].h
                    update_barrier!(model.barrier, h_max_tmp)
                else
                    update_barrier!(model.barrier, model.barrier.elements[model.state.Uk_frame_center].h)
                end

            end
            #  update_barrier!(model.barrier, model.barrier.elements[model.state.Uk_frame_center].h)
        end

        model.stats.nno_successes += 1
    elseif model.state.last_success == PARTIAL_SUCCESS

        # Update the barrier threshold. This follows the implementation of Nomad 3. The threshold h_max
        # is not updated according to the cache but only according to the set of infeasible non dominated points Uᵏ.
        if !isnothing(model.state.h_max)
            below_hxi_elements = filter(elt -> elt.h < model.barrier.elements[model.state.Uk_frame_center].h,
                                        model.barrier.elements[model.barrier.within_Uk .== true])
            if !isempty(below_hxi_elements)
                update_barrier!(model.barrier, maximum(map(elt -> elt.h, below_hxi_elements)))
            end
        end

        model.stats.npartial_successes += 1
    else # Full success

        # update the barrier threshold
        if !isnothing(model.state.h_max)
            below_hmax_elements = filter(elt -> elt.h < model.state.h_max,
                                         model.barrier.elements[model.barrier.within_Uk .== true])
            if !isempty(below_hmax_elements)
                h_max_tmp = maximum(map(elt-> elt.h, below_hmax_elements))
                if h_max_tmp > model.barrier.elements[model.state.Uk_frame_center].h
                    update_barrier!(model.barrier, h_max_tmp)
                else
                    update_barrier!(model.barrier, model.barrier.elements[model.state.Uk_frame_center].h)
                end
            end
        end

        # In the case of phase one, update the mesh; otherwise, it was already set before
        if model.state.is_phase_one
            new_incumbent_index = argmin(map(elt -> elt.h, model.barrier.elements[1:model.barrier.last_index]))
            x_parent = model.cache.Vᵏ[model.barrier.parent_indexes[new_incumbent_index]].inputs
            enlarge_frame_size!(model.barrier.meshes[new_incumbent_index],
                                model.cache.Vᵏ[new_incumbent_index].inputs - x_parent)
        end

        model.stats.nfull_successes += 1
    end

    @assert model.barrier.last_index == model.cache.last_index

    # Reset state
    model.state.last_success = NO_SUCCESS
end

function solve!(model::MadsModel, x0::Vector{Vector{Float64}})
    init_phase!(model, x0)

    init_strategy = begin
        if model.options.use_penalty_approach
            "Use penalty approach."
        elseif model.state.is_phase_one
            "No available starting points: start phase one."
        else
            "Use filter approach."
        end
    end

    model.options.display && println("bb_calls | iteration | id_poll_center | nb of Fk | nb of Ik  |  success type  | hmax")
    model.options.display && println("         |           |                | elements | elements  |                |     ")
    model.options.display && println("-----------------------------------------------------------------------------------")
    model.options.display && println(init_strategy)

    while model.state.stop_reason in [UNDEFINED, STOP_IF_FEASIBLE]
        model.stats.niterations += 1
        set_frame_centers_and_hvalues!(model)
        search_and_poll!(model)
        if model.state.stop_reason == STOP_IF_FEASIBLE
            model.options.display && println("Feasible point found: end of phase one.")
        end
        last_success = model.state.last_success
        update!(model)
        model.options.display && @printf("%5d %10d %10d,%-6d %10d %10d %22s %.9f\n",
                                         model.stats.neval_bb,
                                         model.stats.niterations,
                                         #  model.state.Fk_frame_center,
                                         #  model.state.Uk_frame_center,
                                         model.state.ordered_frame_centers[1],
                                         model.state.ordered_frame_centers[2],
                                         length(get_Fk(model.barrier)),
                                         length(get_Ik(model.barrier)),
                                         last_success,
                                         model.barrier.h_max)
    end
    pb_name = isempty(model.bbproblem.meta.name) ? "Not available" : model.bbproblem.meta.name
    model.options.display && println("-----------------------------------------------------------------------------------")
    println("Summary:")
    println("pb name= ", pb_name)
    println("stop reason= ", model.state.stop_reason)
    println("bb_evals= ", model.stats.neval_bb)
    println("nb_iterations= ", model.stats.niterations)
    println("nb_full_successes= ", model.stats.nfull_successes)
    println("nb_partial_successes= ", model.stats.npartial_successes)
    println("nb_failed_iterations= ", model.stats.nno_successes)
    println("nb_opportunistic_iterations= ", model.stats.nopportunistic_triggers)
    println("nb_cache_hits= ", model.stats.ncache_hits)
    println("nb feasible pts founds= ", length(get_Fk(model.barrier)))
    println("nb infeasible pts founds= ", length(get_Ik(model.barrier)))
    return model.state.stop_reason
end
