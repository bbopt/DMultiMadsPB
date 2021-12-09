export Barrier, insertion_status

export add_feasible!, add_infeasible!, extent, frame_centers, get_Fk, get_Ik, get_Uk, update_barrier!, save_pf_values

const insertion_status = [:dominates, :extends, :improves, :no_improvement, :is_dominated]

mutable struct Barrier

    dims::Tuple{Int, Int} # dimensions of meshes + dimensions of vectors
    h_max::Float64 # h-barrier value

    elements::Vector{OVector}
    meshes::Vector{GranularMesh}
    parent_indexes::Vector{Int}

    within_Fk :: Vector{Bool} # Set of best feasible points
    within_Uk :: Vector{Bool} # Set of current infeasible non dominated points
    within_Ik :: Vector{Bool} # Set of best infeasible non dominated points

    max_size :: Int # Maximum number of objective vectors to save
    last_index :: Int

    function Barrier(dims::Tuple{Int,Int}, max_size::Int; h_max::Float64 = Inf)
        if (dims[1] <= 0) || (dims[2] <= 0)
            error("Non coherent dimensions")
        end
        if h_max < 0
            error("h_max must be positive")
        end
        if max_size <= 0
            error("The barrier cannot contain a negative number of elements")
        end
        return new(dims, h_max,
                   Array{OVector}(undef, max_size),
                   Array{GranularMesh}(undef, max_size),
                   zeros(max_size),
                   falses(max_size), falses(max_size), falses(max_size),
                   max_size, 0)
    end

end

function check_dimension(b::Barrier, v::OVector, m::GranularMesh)
    if b.dims[1] != m.n
        error("Granular mesh and barrier do not have compatible dimensions")
    end
    if b.dims[2] != length(v.f)
        error("Objective vector and barrier do not have compatible dimensions")
    end
end

function check_is_full(b::Barrier)
    if b.last_index == b.max_size
        error("Maximum size of barrier reached: cannot add element")
    end
end

function get_Fk(b::Barrier)
    return @view(b.elements[b.within_Fk .== true])
end

function get_Ik(b::Barrier)
    return @view(b.elements[b.within_Ik .== true])
end

function get_Uk(b::Barrier)
    return @view(b.elements[b.within_Uk .== true])
end

function add_feasible!(b:: Barrier, v::OVector, m::GranularMesh)::Symbol

    # Preliminary checks
    check_dimension(b, v, m)
    check_is_full(b)

    if !isFeasible(v)
        error("Trying to insert an infeasible element into the set of feasible points.")
    end

    Fk_current_elements = @view(b.elements[b.within_Fk .== true])
    nb_elements_in_Fk = length(Fk_current_elements)

    # Empty set
    if nb_elements_in_Fk == 0
        b.last_index += 1
        b.elements[b.last_index] = v
        b.meshes[b.last_index] = deepcopy(m)
        b.within_Fk[b.last_index] = true
        return :improves
    end

    insert = true
    insertion_flag = :improves

    # Check if v dominates an element of Fk
    for index in 1:nb_elements_in_Fk
        corresponding_index = parentindices(Fk_current_elements)[1][index]
        comp_flag = compare(v, b.elements[corresponding_index])
        if comp_flag == :dominating
            b.within_Fk[corresponding_index] = false
            insertion_flag = :dominates
        elseif comp_flag in [:dominated, :equal]
            insertion_flag = :is_dominated
            insert = false
            break
        end
    end

    # Check extension
    if extends(b, v, m, true)
        insertion_flag = :extends
    end

    # Add the point into the barrier ...
    b.last_index += 1
    b.elements[b.last_index] = v
    b.meshes[b.last_index] = deepcopy(m)

    # ... and into the set of feasible points
    if insert
        b.within_Fk[b.last_index] = true
    end

    return insertion_flag
end

function add_infeasible!(b::Barrier, v::OVector, m::GranularMesh)::Symbol

    # Preliminary check
    check_dimension(b, v, m)
    check_is_full(b)

    if isFeasible(v)
        error("Trying to insert a feasible point into the set of infeasible points")
    end

    # Reject the point if above the threshold.
    # Note that this last one is still inserted into the barrier.
    if v.h >= b.h_max
        b.last_index += 1
        b.elements[b.last_index] = v
        b.meshes[b.last_index] = deepcopy(m)
        return :is_rejected
    end

    Uk_current_elements = @view(b.elements[b.within_Uk .== true])
    nb_elements_in_Uk = length(Uk_current_elements)
    prev_nb_best_inf_pts = count(b.within_Ik[1:b.last_index] .== true)

    # Empty set
    if nb_elements_in_Uk == 0
        b.last_index += 1
        b.elements[b.last_index] = v
        b.meshes[b.last_index] = deepcopy(m)
        b.within_Uk[b.last_index] = true
        b.within_Ik[b.last_index] = true
        return :improves
    end

    insert = true

    # Check if v dominates an element of Uk
    for index in 1:nb_elements_in_Uk
        corresponding_index = parentindices(Uk_current_elements)[1][index]
        comp_flag = compare(v, b.elements[corresponding_index])
        if comp_flag == :dominating
            b.within_Uk[corresponding_index] = false
            b.within_Ik[corresponding_index] = false
        elseif comp_flag in [:dominated, :equal]
            insert = false
            break
        end
    end

    insertion_flag = :improves

    # Check extension
    if extends(b, v, m, false)
        insertion_flag = :extends
    end

    # Add the point into the barrier ...
    b.last_index += 1
    b.elements[b.last_index] = v
    b.meshes[b.last_index] = deepcopy(m)

    # ... and into the set of filter points
    if insert
        b.within_Uk[b.last_index] = true
    else
        return :is_dominated
    end

    insertion_best_flag = :is_dominated

    # Update set of best infeasible points
    if insert
        insertion_best_flag = _update_Ik_set_after_insertion!(b)
    end

    # Set corresponding flags according to the situation
    if insertion_flag == :extends
        return insertion_flag
    end

    if insertion_best_flag == :is_dominated
        # At this moment, the point has been successfully inserted into the filter
        return :no_improvement
    else
        # The number of best infeasible points may have been already reduced
        # by the insertion into the filter, hence this check.
        if count(b.within_Ik) <= prev_nb_best_inf_pts
            return :dominates
        else
            return insertion_best_flag
        end
    end
end

function _private_compare_Ik_elements(v1::OVector, v2::OVector)::Symbol
    isbetter = false
    isworse = false
    for (f1, f2) in zip(v1.f, v2.f)
        if (f1 < f2)
            isbetter = true
        end
        if (f2 < f1)
            isworse = true
        end
        if (isworse && isbetter)
            break
        end
    end
    if isworse
        if isbetter
            return :nondominated
        else
            return :dominated
        end
    else
        if isbetter
            return :dominating
        else
            return :equal
        end
    end
end

# This function must be called when the last point has been inserted.
function _update_Ik_set_after_insertion!(b::Barrier)::Symbol
    candidate = b.elements[b.last_index]

    domination_flags = falses(b.last_index)

    insert = true
    insertion_flag = :improves

    # Check the candidate dominates an element in Ik
    for (index, is_in_Ik) in enumerate(@view(b.within_Ik[1:b.last_index-1]))
        if is_in_Ik
            comp_flag = _private_compare_Ik_elements(candidate, b.elements[index])
            if comp_flag == :dominating
                domination_flags[index] = true
                insertion_flag = :dominates
            elseif comp_flag in [:dominated, :equal]
                insertion_flag = :is_dominated
                insert = false
                break
            end
        end
    end

    # Update the set of best infeasible points
    Ik_view = @view(b.within_Ik[1:b.last_index])
    Ik_view[domination_flags .== true] .= false
    Ik_view[b.last_index] = insert

    return insertion_flag
end

function extends(b :: Barrier, v :: OVector, m::GranularMesh, is_feasible::Bool)::Bool
    check_dimension(b, v, m)
    tmp_best_pts = begin
        if is_feasible
            @view b.elements[b.within_Fk]
        else
            @view b.elements[b.within_Uk]
        end
    end
    if isempty(tmp_best_pts)
        return false
    end

    # Ideal vector of the best points
    ideal_v = Inf * ones(b.dims[2], 1)
    for elt in tmp_best_pts
        ideal_v = min.(elt.f, ideal_v)
    end
    return any(v.f .< ideal_v)
end

function update_barrier!(b ::Barrier,  h_max:: Float64)
    if h_max < 0
        error("h_max cannot be negative")
    end
    b.h_max = h_max

    filter_flags = falses(b.last_index)

    for (index, is_in_Uk) in enumerate(b.within_Uk[1:b.last_index])
        if is_in_Uk
            if b.elements[index].h > b.h_max
                filter_flags[index] = true
            end
        end
    end

    # Remove all Uk elements above the threshold
    Uk_view = @view(b.within_Uk[1:b.last_index])
    Uk_view[filter_flags .== true] .= false

    # Remove all Ik elements above the threshold
    Ik_view = @view(b.within_Ik[1:b.last_index])
    Ik_view[filter_flags .== true] .= false

    # Reinsert potential new non dominated points into the set Ik
    _update_Ik_after_h_max_setting!(b)
end

# This function has to be called after setting h_max.
function _update_Ik_after_h_max_setting!(b::Barrier)
    for (index, is_in_Uk) in enumerate(b.within_Uk[1:b.last_index])
        if is_in_Uk && !b.within_Ik[index]

            # Check if element indexed by index can be inserted in the set Ik
            insert = true
            for (index_2, is_in_Uk_2) in enumerate(@view(b.within_Uk[1:b.last_index]))
                if (is_in_Uk_2) && (index != index_2)
                    comp_flag = _private_compare_Ik_elements(b.elements[index], b.elements[index_2])
                    if comp_flag == :dominating
                        b.within_Ik[index_2] = false
                    elseif comp_flag in [:dominated, :equal]
                        insert = false
                        break
                    end
                end
            end
            b.within_Ik[index] = insert
        end
    end
end

# Return the indices of the frame centers corresponding to the barrier `elements` attribute
function frame_centers(b::Barrier, w::Int; use_doM_selection::Bool=true)

    feasible_index = feasible_frame_center(b, w)
    infeasible_index = begin
        if iszero(feasible_index)
            infeasible_frame_center(b)
        else
            # TODO : change the distance selection
            function dom_distance(z1::Vector{Float64}, z2::Vector{Float64})::Float64
                if any(z1 .< z2)
                    return norm(z2 - z1)
                else
                    return -norm(z2 - z1)
                end
            end

            if use_doM_selection
                # Get the Ik point below the h_max feasible threshold which has
                # maximal dominance move and has not been assigned yet.
                tmp_ind = 0
                tmp_dist = -Inf
                for ind in 1:b.last_index
                    if b.within_Ik[ind]
                        v = b.elements[ind]
                        tmp_dom_distance = minimum(elt -> sum(elt.f - min.(v.f, elt.f)), get_Fk(b))
                        if v.h <= b.h_max && tmp_dom_distance > tmp_dist
                            tmp_ind = ind
                            tmp_dist = tmp_dom_distance
                        end
                    end
                end

                # all Ik points are dominated by an element of Fk
                if tmp_dist == 0
                    # In this case, get the Ik point below the h_max feasible threshold which has
                    # minimal dominance move, when considered a maximization problem, and has not
                    # been assigned yet.
                    tmp_ind =  0
                    tmp_dist = Inf
                    for ind in 1:b.last_index
                        if b.within_Ik[ind]
                            v = b.elements[ind]
                            # TODO to check this distance
                            tmp_dom_distance = minimum(elt -> sum(v.f - min.(v.f, elt.f)), get_Fk(b))
                            if v.h <= b.h_max && tmp_dom_distance < tmp_dist
                                tmp_ind = ind
                                tmp_dist = tmp_dom_distance
                            end
                        end
                    end
                end
                tmp_ind
            else
                # Get the Ik point below the h_max feasible threshold which has minimal value.
                tmp_ind = 0
                tmp_dist = Inf
                for ind in 1:b.last_index
                    if b.within_Ik[ind]
                        v = b.elements[ind]
                        tmp_dom_distance = dom_distance(b.elements[feasible_index].f, v.f)
                        if v.h <= b.h_max && tmp_dom_distance < tmp_dist
                            tmp_ind = ind
                            tmp_dist = tmp_dom_distance
                        end
                    end
                end
                tmp_ind
            end
        end
    end

    return (feasible = feasible_index,
            infeasible = infeasible_index)
end

function feasible_frame_center(b::Barrier, w::Int)::Int

    if isempty(get_Fk(b))
        return 0
    end

    Fk_indexes = parentindices(get_Fk(b))[1]

    # Get maximum mesh size value
    Fk_Δ_max = 0.0
    for ind in Fk_indexes
        Δ_val = norm(get_frame_size_parameter(b.meshes[ind]), Inf)
        Fk_Δ_max = max(Fk_Δ_max, Δ_val)
    end

    # Select candidates
    Fk_selected_indexes = Int[]
    for ind in Fk_indexes
        Fk_elt_mesh = b.meshes[ind]
        Δ_val = norm(get_frame_size_parameter(Fk_elt_mesh), Inf)
        if (10.0^(-w) * Fk_Δ_max) <= Δ_val && !check_mesh_for_stopping(Fk_elt_mesh)
            push!(Fk_selected_indexes, ind)
        end
    end

    # The selection must always work
    if isempty(Fk_selected_indexes)
        return Fk_indexes[1]
    end

    # Only one point
    if length(Fk_selected_indexes) == 1
        return Fk_selected_indexes[1]
    end

    # Two points in the barrier
    if length(Fk_selected_indexes) == 2 && length(Fk_indexes) == 2
        v1f = b.elements[Fk_selected_indexes[1]].f
        v2f = b.elements[Fk_selected_indexes[2]].f
        if norm(v1f, Inf) > norm(v2f, Inf)
            return Fk_selected_indexes[1]
        else
            return Fk_selected_indexes[2]
        end
    end

    # More than two points
    frame_ind = 0
    maximum_gap = -1.0 # negative to deal with the case where two points in Fk_selected_indexes

    for obj in 1: b.dims[2]
        fvalues = [(b.elements[ind].f[obj], ind) for ind in Fk_indexes]
        fvalues = sort(fvalues)

        # Get extreme points according to one objective
        fmin = fvalues[1][1]
        fmax = fvalues[end][1]

        # Can happen for exemple when we have several minima or for more than three objectives
        if fmin == fmax
            fmin = 0.0
            fmax = 1.0
        end

        # Intermediate points
        for i in 2:length(fvalues)-1
            current_gap = (fvalues[i+1][1] - fvalues[i-1][1]) / (fmax - fmin)
            if fvalues[i][2] in Fk_selected_indexes && current_gap >= maximum_gap
                maximum_gap = current_gap
                frame_ind = fvalues[i][2]
            end
        end

        # Extreme points
        current_gap = 2 * (fvalues[end][1] - fvalues[end-1][1]) / (fmax - fmin)
        if fvalues[end][2] in Fk_selected_indexes && current_gap > maximum_gap
            maximum_gap = current_gap
            frame_ind = fvalues[end][2]
        end

        current_gap = 2 * (fvalues[2][1] - fvalues[1][1]) / (fmax - fmin)
        if fvalues[1][2] in Fk_selected_indexes && current_gap > maximum_gap
            maximum_gap = current_gap
            frame_ind = fvalues[1][2]
        end
    end

    return frame_ind
end

function infeasible_frame_center(b::Barrier)::Int

    if isempty(get_Uk(b))
        return 0
    end

    Ik_indexes = parentindices(get_Ik(b))[1]
    Ik_Δ_min = let xI_ind = argmin(map(elt -> elt.h, b.elements[Ik_indexes]))
        norm(get_frame_size_parameter(b.meshes[Ik_indexes[xI_ind]]), Inf)
    end

    Ik_selected_indexes = Int[]
    for ind in Ik_indexes
        Δ_val = norm(get_frame_size_parameter(b.meshes[ind]), Inf)
        if Ik_Δ_min <= Δ_val
            push!(Ik_selected_indexes, ind)
        end
    end

    if isempty(Ik_selected_indexes)
        return 0
    end

    # Only one point
    if length(Ik_selected_indexes) == 1
        return Ik_selected_indexes[1]
    end

    # Two points into the barrier
    if length(Ik_selected_indexes) == 2 && length(Ik_indexes) == 2
        v1f = b.elements[Ik_selected_indexes[1]].f
        v2f = b.elements[Ik_selected_indexes[2]].f
        if norm(v1f, Inf) > norm(v2f, Inf)
            return Ik_selected_indexes[1]
        else
            return Ik_selected_indexes[2]
        end
    end

    # More than two points
    frame_ind = 0
    maximum_gap = -1.0 # negative to deal with the case where two points in Ik_selected_indexes

    for obj in 1: b.dims[2]
        fvalues = [(b.elements[ind].f[obj], ind) for ind in Ik_indexes]
        fvalues = sort(fvalues)

        # Get extreme points according to one objective
        fmin = fvalues[1][1]
        fmax = fvalues[end][1]

        # Can happen for exemple when we have several minima or for more than three objectives
        if fmin == fmax
            fmin = 0.0
            fmax = 1.0
        end

        # Intermediate points
        for i in 2:length(fvalues)-1
            current_gap = (fvalues[i+1][1] - fvalues[i-1][1]) / (fmax - fmin)
            if fvalues[i][2] in Ik_selected_indexes && current_gap >= maximum_gap
                maximum_gap = current_gap
                frame_ind = fvalues[i][2]
            end
        end

        # Extreme points
        current_gap = 2 * (fvalues[end][1] - fvalues[end-1][1]) / (fmax - fmin)
        if fvalues[end][2] in Ik_selected_indexes && current_gap > maximum_gap
            maximum_gap = current_gap
            frame_ind = fvalues[end][2]
        end

        current_gap = 2 * (fvalues[2][1] - fvalues[1][1]) / (fmax - fmin)
        if fvalues[1][2] in Ik_selected_indexes && current_gap > maximum_gap
            maximum_gap = current_gap
            frame_ind = fvalues[1][2]
        end
    end

    return frame_ind
end

# Return the extent of the Pareto front (a customized one)
function extent(b::Barrier)::Float64
    if isempty(get_Fk(b))
        return 0
    end

    Fk_indexes = parentindices(get_Fk(b))[1]

    extent_val = 0

    for obj in 1: b.dims[2]
        fvalues = [b.elements[ind].f[obj] for ind in Fk_indexes]
        fvalues = sort(fvalues)

        # Get extreme points according to one objective
        fmin = fvalues[1]
        fmax = fvalues[end]

        # Can happen for exemple when we have several minima or for more than three objectives
        if fmin == fmax
            extent_val += abs(fmin)
        else
            extent_val += abs(fmax - fmin)
        end
    end
    return extent_val
end

function save_pf_values(b :: Barrier, filename :: String)
    open(filename, "w") do io
        # write values of Pareto front
        for elt in get_Fk(b)
            writedlm(io, transpose(elt.f))
        end
    end
end
