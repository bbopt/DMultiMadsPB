export EvalPoint, Cache, isincache, add_cache!, get_Vk, save_cache

# cache structure and eval point for the problem
struct EvalPoint

    inputs :: Vector{Float64} # Inputs of the black box
    outputs :: Vector{Float64} # Outputs of the black box

    function EvalPoint(inputs, outputs)
        if (length(inputs) < 1) || (length(outputs) < 1)
            error("Non sensical dimensions")
        end

        new(inputs, outputs)
    end
end

# print meta information EvalPoint
import Base.print, Base.show, Base.println
function show(io :: IO, x :: EvalPoint)
    s = @sprintf("Inputs = %s, ", string(x.inputs))
    s *= @sprintf("Outputs = %s\n", string(x.outputs))
    print(io, s)
end

function print(io :: IO, x :: EvalPoint)
    @printf(io, "Inputs = %s, ", string(x.inputs))
    @printf(io, "Outputs = %s\n", string(x.outputs))
end

# this structure keeps all non-redundant points evaluated during the iterations
mutable struct Cache

    dims:: Tuple{Int, Int} # precise the dimensions of the inputs and outputs it can take
    Vᵏ  :: Vector{EvalPoint} # storing of the data

    max_size :: Int
    last_index :: Int

    function Cache(dimensions::Tuple{Int, Int}, max_size::Int)
        if dimensions[1] <= 0 || dimensions[2] <= 0
            error("Non sensical dimensions")
        end
        if max_size <= 0
            error("A cache cannot have a negative number of elements")
        end
        new(dimensions, Array{EvalPoint}(undef, max_size), max_size, 0)
    end

end

# return true if the point is in cache, false otherwise
function isincache(c :: Cache, x :: Vector{Float64})::Bool
    if length(x) != c.dims[1]
        return false
    end

    in_cache = false
    for elt in @view c.Vᵏ[1:c.last_index]
        if maximum(abs.(elt.inputs - x)) <= 10^(-9)
            in_cache = true
            break
        end
    end
    return in_cache
end

# return true if the add_cache function works, false otherwise
function add_cache!(c :: Cache, x :: EvalPoint)::Bool
    if (c.max_size == c.last_index)
        error("The cache is full: cannot add element")
    end
    isOk = isincache(c, x.inputs)
    if isOk == false
        if (length(x.inputs) == c.dims[1]) && (length(x.outputs) == c.dims[2])
            c.last_index += 1
            c.Vᵏ[c.last_index] = x
            return true
        end
    end
    return false
end

function get_Vk(c::Cache)
    return @view c.Vᵏ[1:c.last_index]
end

function save_cache(c :: Cache, filename :: String)
    open(filename, "w") do io
        # write dimensions of the cache
        writedlm(io, transpose(collect(c.dims)))
        # write inputs and outputs in the file
        for elt in @view c.Vᵏ[1:c.last_index]
            writedlm(io, transpose([elt.inputs; elt.outputs]))
        end
    end
end

# print meta information on cache
function print(io :: IO , c :: Cache)
    @printf(io, "Dimensions elements = %s\n", string(c.dims))
    @printf(io, "Number elements = %d\n\n", length(c.Vᵏ))
    for elt in @view c.Vᵏ[1:c.last_index]
        println(io, elt)
    end
end
