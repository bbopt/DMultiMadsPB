export AbstractStepSet, AbstractPollSet, AbstractSearchSet,
PollSet1, PollSet2, PollSet2N, PollSetNp1,
SpeculativeSearchSet
#, IterationCenterSet

export generate_candidates, generate_ordered_candidates, generate_ordered_candidates_per_step

# generic type for all step containers
abstract type AbstractStepSet end

# generic type for all poll steps
abstract type AbstractPollSet <: AbstractStepSet end

# generic type for all search steps
abstract type AbstractSearchSet <: AbstractStepSet end

# 1- Poll sets

# classical poll with 2n Orthomads directions
struct PollSet2N <: AbstractPollSet
    mesh :: GranularMesh
    poll_center :: Vector{Float64}
    last_success_direction:: Union{Nothing, Vector{Float64}}

    function PollSet2N(m :: GranularMesh, poll_center::Vector{Float64}; 
                       last_success_direction::Union{Nothing, Vector{Float64}} = nothing)
        @lencheck m.n poll_center
        if last_success_direction != nothing
            @lencheck m.n last_success_direction
        end

        # TODO make deep copy of mesh?
        new(m, poll_center, last_success_direction)
    end

end

function generate_candidates(p :: PollSet2N, rng::AbstractRNG)::Array{Float64, 2}
    # set direction on the unit sphere
    dir = randn(rng, (p.mesh.n))
    dir /= norm(dir)

    # householder transformation
    H = I - 2 * dir * transpose(dir)

    # scale and project directions on the mesh
    for j in 1:size(H)[2]
        H[:, j] = scale_and_project_on_mesh(p.mesh, H[:, j])
    end

    return p.poll_center .+ [H -H]

end

# reduced poll with n+1 Orthomads directions
struct PollSetNp1 <: AbstractPollSet
    mesh :: GranularMesh
    poll_center :: Vector{Float64}
    last_success_direction:: Union{Nothing, Vector{Float64}}

    function PollSetNp1(m :: GranularMesh, poll_center::Vector{Float64}; 
                        last_success_direction::Union{Nothing, Vector{Float64}} = nothing)
        @lencheck m.n poll_center
        if last_success_direction != nothing
            @lencheck m.n last_success_direction
        end

        # TODO make deep copy of mesh?
        new(m, poll_center, last_success_direction)
    end

end

function generate_candidates(p :: PollSetNp1, rng::AbstractRNG)::Array{Float64, 2}
    # set direction on the unit sphere
    dir = randn(rng, (p.mesh.n))
    dir /= norm(dir)

    # householder transformation
    H = I - 2 * dir * transpose(dir)

    # scaling and project directions on the mesh
    for j in 1:size(H)[2]
        H[:, j] = scale_and_project_on_mesh(p.mesh, H[:, j])
    end

    # succ-neg strategy as explained in the paper for orthomads
    if p.last_success_direction != nothing
        H_succ_neg = copy(H)
        for j in 1:size(H)[2]
            if dot(p.last_success_direction, H[:, j]) >= 0
                H_succ_neg[:, j] = H[:, j]
            else
                H_succ_neg[:, j] = -H[:, j]
            end
        end
        return p.poll_center .+ [H_succ_neg -sum(H_succ_neg, dims=2)]
        # 2 n success strategy
    else
        return p.poll_center .+ [H -H]
    end

    # scale and project directions on the mesh
    return p.poll_center .+ [H -H]

end

# poll with 2 directions
struct PollSet2 <: AbstractPollSet
    mesh :: GranularMesh
    poll_center :: Vector{Float64}
    last_success_direction:: Union{Nothing, Vector{Float64}}

    function PollSet2(m :: GranularMesh, poll_center::Vector{Float64};
                      last_success_direction::Union{Nothing, Vector{Float64}} = nothing)
        @lencheck m.n poll_center
        if last_success_direction != nothing
            @lencheck m.n last_success_direction
        end

        # TODO make deep copy of mesh?
        new(m, poll_center, last_success_direction)
    end

end

function generate_candidates(p :: PollSet2, rng::AbstractRNG)::Array{Float64, 2}
    # set direction on the unit sphere
    dir = randn(rng, (p.mesh.n))
    dir /= norm(dir)

    dirs = [dir -dir]

    # scale and project directions on the mesh
    for j in 1:size(dirs)[2]
        dirs[:, j] = scale_and_project_on_mesh(p.mesh, dirs[:, j])
    end

    return p.poll_center .+ dirs
end

# poll with 1 direction
struct PollSet1 <: AbstractPollSet
    mesh :: GranularMesh
    poll_center :: Vector{Float64}
    last_success_direction :: Union{Nothing, Vector{Float64}}

    function PollSet1(m :: GranularMesh, poll_center::Vector{Float64};
                      last_success_direction::Union{Nothing, Vector{Float64}} = nothing)
        @lencheck m.n poll_center
        if last_success_direction != nothing
            @lencheck m.n last_success_direction
        end

        # TODO make deep copy of mesh?
        new(m, poll_center, last_success_direction)
    end

end

function generate_candidates(p :: PollSet1, rng::AbstractRNG)::Array{Float64, 2}
    # set direction on the unit sphere
    dir = randn(rng, (p.mesh.n))
    dir /= norm(dir)

    # scale and project directions on the mesh
    dir = scale_and_project_on_mesh(p.mesh, dir)

    return reshape(p.poll_center .+ dir, length(p.poll_center), 1)
end

# generate ordered candidates by last success direction
function generate_ordered_candidates(p :: AbstractPollSet, rng::AbstractRNG)::Array{Float64, 2}
    candidates = generate_candidates(p, rng)
    if p.last_success_direction != nothing

        function get_angle(dir1, dir2)
            val = dot(dir1, dir2) / (norm(dir1) * norm(dir2))
            if abs(val) > 1
                return acos(1)
            else
                return acos(val)
            end
        end
        candidates = sortslices(candidates, dims=2,
                                by=elt->get_angle(elt - p.poll_center, p.last_success_direction))
    end
    return candidates
end

# 2- Search sets

# speculative search step
struct SpeculativeSearchSet <: AbstractSearchSet
    mesh :: GranularMesh
    mesh_center :: Vector{Float64}
    last_success_direction:: Union{Nothing, Vector{Float64}}

    function SpeculativeSearchSet(m :: GranularMesh, mesh_center :: Vector{Float64};
                                  last_success_direction:: Union{Nothing, Vector{Float64}}=nothing)
        @lencheck m.n mesh_center
        if last_success_direction != nothing
            @lencheck m.n last_success_direction 
        end

        # TODO make deep copy of mesh?
        new(m, mesh_center, last_success_direction)

    end

end

function generate_candidates(s::SpeculativeSearchSet, rng::AbstractRNG)
    if s.last_success_direction == nothing
        return Array{Float64}(undef, s.mesh.n, 0)
    else
        factor = Inf
        for i in 1:s.mesh.n
            if s.last_success_direction[i] != 0
                factor = min(factor, get_frame_size_parameter(s.mesh)[i] / abs(s.last_success_direction[i]))
            end
        end
        candidates = s.mesh_center .+ factor * s.last_success_direction
        return reshape(candidates, length(candidates), 1)
    end
end

function generate_ordered_candidates(s::SpeculativeSearchSet, rng::AbstractRNG)
    return generate_candidates(s, rng)
end
