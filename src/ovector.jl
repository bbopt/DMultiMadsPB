export OVector, ≦, <=, <, compare, isEqual, isFeasible

# Outcome Vector
struct OVector

    f :: Vector{Float64} # objective function
    h :: Float64 # h value

    function OVector(f, h)
        if isempty(f)
            error("objective function cannot be empty")
        end
        if h < 0
            error("h value cannot be negative")
        end
        new(f, h)
    end

end

# print meta information OVector
import Base.print, Base.show, Base.println
function show(io :: IO, v :: OVector)
    s = @sprintf("f = %s, ", string(v.f))
    s *= @sprintf("h = %s, ", string(v.h))
    print(io, s)
end

function print(io :: IO, v :: OVector)
    @printf(io, "f = %s, ", string(v.f))
    @printf(io, "h = %s\n", string(v.h))
end

# check which category does it belong
function isFeasible(v :: OVector)::Bool
    if v.h == 0
        return true
    else
        return false
    end
end

# return dominance comparison in at least m scalar comparisons
function compare(v1::OVector, v2::OVector)::Symbol
    if length(v1.f) != length(v2.f) || !((isFeasible(v1) && isFeasible(v2)) || (!isFeasible(v1) && !isFeasible(v2)))
        return :undefined
    end

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
    if !(isworse && isbetter)
        if (v1.h < v2.h)
            isbetter = true
        end
        if (v2.h < v1.h)
            isworse = true
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

# relation orders

# weak dominance
function ≦(v1 :: OVector, v2 :: OVector)::Bool
    return compare(v1, v2) in [:equal, :dominating]
end

# strict dominance
import Base.<
function <(v1 :: OVector, v2 :: OVector)::Bool
    if length(v1.f) != length(v2.f)
        return false
    end
    strict_dom = false;
    if isFeasible(v1) && isFeasible(v2)
        strict_dom = all(v1.f .< v2.f)
    elseif !isFeasible(v1) && !isFeasible(v2)
        strict_dom = all(v1.f .< v2.f) && (v1.h < v2.h)
    end
    return strict_dom
end

# dominance
import Base.<=
function <=(v1::OVector, v2::OVector)::Bool
    return compare(v1, v2) == :dominating
end

# equality
function isEqual(v1::OVector, v2::OVector)
    return compare(v1, v2) == :equal
end

# relation orders :
# by convention, a feasible vector does not dominates an infeasible vector
# weak dominance
#  function ≦(v1 :: OVector, v2 :: OVector)::Bool
#      if length(v1.f) != length(v2.f)
#          return false
#      end
#      weakly_dom = true && (isFeasible(v1) && isFeasible(v2)) || (!isFeasible(v1) && !isFeasible(v2))
#      if isFeasible(v1) && isFeasible(v2)
#          for (f1, f2) in zip(v1.f, v2.f)
#              if f1 > f2
#                  weakly_dom = false
#                  break
#              end
#          end
#          #  weakly_dom = all(v1.f .<= v2.f)
#      elseif !isFeasible(v1) && !isFeasible(v2)
#          for (f1, f2) in zip(v1.f, v2.f)
#              if f1 > f2
#                  weakly_dom = false
#                  break
#              end
#          end
#          weakly_dom = weakly_dom && (v1.h <= v2.h)
#      end
#      return weakly_dom
#  end
#
#  import Base.<
#  # strict dominance
#  function <(v1 :: OVector, v2 :: OVector)::Bool
#      if length(v1.f) != length(v2.f)
#          return false
#      end
#      strict_dom = false;
#      if isFeasible(v1) && isFeasible(v2)
#          strict_dom = all(v1.f .< v2.f)
#      elseif !isFeasible(v1) && !isFeasible(v2)
#          strict_dom = all(v1.f .< v2.f) && (v1.h < v2.h)
#      end
#      return strict_dom
#  end
#
#  # dominance
#  import Base.<=
#  function <=(v1 :: OVector, v2 :: OVector)::Bool
#      if length(v1.f) != length(v2.f)
#          return false
#      end
#      isbetter = false
#      isworse = false
#
#      if isFeasible(v1) && isFeasible(v2)
#          for (f1, f2) in zip(v1.f, v2.f)
#              if (f1 < f2)
#                  isbetter = true
#              end
#              if (f2 < f1)
#                  isworse = true
#              end
#              if (isworse && isbetter)
#                  break
#              end
#          end
#          #  dominance = all(v1.f .<= v2.f) && any(v1.f .< v2.f)
#          if isworse
#              return false
#          else
#              if isbetter
#                  return true
#              else
#                  return false
#              end
#          end
#      elseif !isFeasible(v1) && !isFeasible(v2)
#          dominance = (all(v1.f .<= v2.f) && (v1.h < v2.h)) || (all(v1.f .<= v2.f) && any(v1.f .< v2.f) && (v1.h <= v2.h))
#          return dominance
#      end
#      return (isFeasible(v1) && isFeasible(v2)) || (!isFeasible(v1) && !isFeasible(v2))
#  end
#
#  # equality
#  function isEqual(v1 :: OVector, v2:: OVector)::Bool
#      return (v1.f ≈ v2.f) && (v1.h ≈ v2.h)
#  end
