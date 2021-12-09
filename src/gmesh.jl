# Implementation of the granular as exposed in the article
#
# The Mesh Adaptive Direct Search Algorithm for Granular and Discrete Variables
# Charles Audet, Sébastien Le Digabel, and Christophe Tribes
# SIAM Journal on Optimization 2019 29:2, 1164-1189
#
# implementation of the granular mesh : perform better for the mesh adaptive direct search algorithm
# all variables are continuous in our case
# The notations follow the paper

export
    GranularMesh,
    get_mesh_size_parameter,
    get_frame_size_parameter,
    enlarge_frame_size!,
    refine_frame_size!,
    scale_and_project_on_mesh,
    project_on_mesh,
    project_on_mesh_and_snap_to_bounds

struct GranularMesh

    n :: Int # dimension of the mesh
    δ_min :: Float64 # minimum mesh size
    b_init :: Vector{Int} # initial poll size exponent

    bᵏ :: Vector{Int} # current poll size exponent
    aᵏ :: Vector{Int} # current poll size mantisse

    function GranularMesh(n :: Int, Δ_init :: Vector{Float64} ;
                   δ_min :: Float64 = 10^(-13))
        if n < 1
            error("Non sensical dimensions")
        end

        if sum(Δ_init .<= 0) > 0
            error("Initial frame size parameter cannot be negative")
        end

        @lencheck n Δ_init

        exponents = trunc.(Int, log10.(abs.(Δ_init)))

        function round_poll_size_mant(mant :: Float64)::Int
            # Round input mant to 1, 2 or 5 for poll_size_mant element
            poll_size_mant = 0
            if mant < 1.5
                poll_size_mant = 1
            elseif ( mant >= 1.5) && (mant < 3.5 )
                poll_size_mant = 2
            else
                poll_size_mant = 5
            end

            return poll_size_mant
        end

        mantisses = [round_poll_size_mant(Δ_init[i] * 10.0^(-exponents[i])) for i in 1:length(Δ_init)]

        new(n, δ_min, copy(exponents), exponents, mantisses)
    end

end

function get_mesh_size_parameter(m :: GranularMesh)::Vector{Float64}
    δ = 10.0 .^ (m.bᵏ - abs.(m.bᵏ - m.b_init))
    # if the mesh size parameter is below the minimum mesh size, return the min mesh size
    return max.(m.δ_min, δ)
end

function get_frame_size_parameter(m :: GranularMesh)::Vector{Float64}
    return m.aᵏ .* 10.0 .^ m.bᵏ
end

function get_ρk(m :: GranularMesh)::Vector{Float64}
    return m.aᵏ .* 10.0 .^ (abs.(m.bᵏ .- m.b_init))
end

function enlarge_frame_size!(m :: GranularMesh, direction::Vector{Float64}; aₜ::Float64 = 0.1 )::Bool
    # a_t is the anisotropy factor
    @lencheck m.n direction

    ischanged = false
    ρ_min = min(Inf, minimum(get_ρk(m)))

    for i in 1:m.n
        if (abs(direction[i]) / get_mesh_size_parameter(m)[i] / get_ρk(m)[i] > aₜ) || ((m.bᵏ[i] < m.b_init[i]) && (get_ρk(m)[i] > ρ_min * ρ_min))

            if m.aᵏ[i] == 1
                m.aᵏ[i] = 2
            elseif m.aᵏ[i] == 2
                m.aᵏ[i] = 5
            else
                m.aᵏ[i] = 1
                m.bᵏ[i] += 1
            end
            ischanged = true
        end
    end
    return true
end

function refine_frame_size!(m :: GranularMesh)

    # compute new poll_size_mant and poll_size_exp
    poll_size_mant = copy(m.aᵏ)
    poll_size_exp = copy(m.bᵏ)

    for i in 1:m.n
        if poll_size_mant[i] == 1
            poll_size_mant[i] = 5
            poll_size_exp[i] -= 1
        elseif poll_size_mant[i] == 2
            poll_size_mant[i] = 1
        else
            poll_size_mant[i] = 2
        end
    end

    # safety check
    δ_old = get_mesh_size_parameter(m)
    for i in 1: m.n
        if m.δ_min < δ_old[i]
            m.aᵏ[i] = poll_size_mant[i]
            m.bᵏ[i] = poll_size_exp[i]
        end
    end

end

function check_mesh_for_stopping(m :: GranularMesh)::Bool
    δ = get_mesh_size_parameter(m)
    if any(m.δ_min .>= δ)
        return true
    else
        return false
    end
end

# used for projecting directions
function scale_and_project_on_mesh(m:: GranularMesh, elt :: Vector{Float64})::Vector{Float64}
    @lencheck m.n elt
    if norm(elt, Inf) == 0
        error("Cannot handle an infinite norm of zeros")
    end
    δ = get_mesh_size_parameter(m)
    rho = get_ρk(m)
    return round.(rho .* elt / norm(elt, Inf)).*δ
end

# the projection is made finding the closest element of the mesh centered at x_center from proj
function project_on_mesh(m :: GranularMesh, proj :: Vector{Float64}, x_center :: Vector{Float64})::Vector{Float64}
    @lencheck m.n proj x_center
    δ = get_mesh_size_parameter(m)
    candidate = zeros(m.n, )
    for i in 1 : m.n
        candidate[i] = δ[i] * round((proj[i] - x_center[i]) / δ[i]) + x_center[i]
    end
    return candidate
end

function project_on_mesh_and_snap_to_bounds(m::GranularMesh,
                                            proj :: Vector{Float64},
                                            x_center :: Vector{Float64},
                                            lb:: Vector{Float64},
                                            ub::Vector{Float64})
    @lencheck m.n lb ub
    if !all(lb .< ub)
        error("Wrong bound constraints")
    end
    if !all(lb .<= x_center) || !all(ub .>= x_center)
        error("mesh center values must satisfy: lb <= x^mesh <= ub")
    end

    # 1- project on the mesh
    candidate = project_on_mesh(m, proj, x_center)

    # 2- snap to bound if necessary
    δ = get_mesh_size_parameter(m)
    snapped_candidate = zeros(m.n, )
    for i in 1: m.n
        if (lb[i] <= candidate[i]) && (candidate[i] <= ub[i])
            snapped_candidate[i] = candidate[i]
        else
            if candidate[i] < lb[i]
                # meshCenter is supposed to be >= lb; normally, this rounding is supposed to be in the box constraints
                snapped_candidate[i] = δ[i] * ceil((lb[i] - x_center[i]) / δ[i]) + x_center[i]
            else
                # ref_value is supposed to be <= ub; normally, this rounding is supposed to be in the box constraints
                snapped_candidate[i] =  δ[i] * floor((ub[i] - x_center[i]) / δ[i]) + x_center[i]
            end
            # as defined in Nomad 3
            if snapped_candidate[i] < lb[i]
                @warn begin
                """
                Warning: snap_to_bounds: Error snapping $(candidate[i]) to lower bound $(lb[i])
                frameCenter = $(x_center[i]), δ = $(δ[i]) : it gave $(snapped_candidate[i]) which is still lower than $(lb[i])
                """
                # TODO Force the snapping ?
                end
            end
            if snapped_candidate[i] > ub[i]
                @warn begin
                """
                Warning: snap_to_bounds: Error snapping $(candidate[i]) to upper bound $(ub[i])
                frameCenter = $(x_center[i]), δ = $(δ[i]) : it gave $(snapped_candidate[i]) which is still higher than $(ub[i])
                """
                # TODO Force the snapping ?
                end
            end
        end
    end

    return snapped_candidate
end

# print meta information gmesh
import Base.print, Base.show, Base.println
function show(io :: IO, m :: GranularMesh)
    s = @sprintf("Dimensions = %d\n", m.n)
    s *= @sprintf("δ_min = %d\n", m.δ_min)
    s *= @sprintf("init_poll_size_exp = %s\n", string(m.b_init))
    s *= @sprintf("poll_size_exp = %s\n", string(m.bᵏ))
    s *= @sprintf("poll_size_mant = %s\n", string(m.aᵏ))
    print(io, s)
end

function print(io :: IO, m :: GranularMesh)
    @printf(io, "Dimensions = %d\n", m.n)
    @printf(io, "δ_min = %d\n", m.δ_min)
    @printf(io, "init_poll_size_exp = %s\n", string(m.b_init))
    @printf(io, "poll_size_exp = %s\n", string(m.bᵏ))
    @printf(io, "poll_size_mant = %s\n", string(m.aᵏ))
    @printf(io, "δᵏ = %s\n", string(get_mesh_size_parameter(m)))
    @printf(io, "Δᵏ = %s\n", string(get_frame_size_parameter(m)))
end
