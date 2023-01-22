module TaylorInversion
using TaylorSeries: Taylor1
using Symbolics: @variables, @acrule, @rule
using Symbolics: Arr, Num, Chain, Fixpoint, Prewalk
using Symbolics: build_function, solve_for, expand, simplify, substitute

struct InverseTaylor{N}
    a::Arr{Num,1}
    A::Arr{Num,1}
    z::Num
    B::Vector{Num}
    function InverseTaylor{N}() where {N}
        @variables z, a[1:N], A[1:N]
        B = zeros(Num, N)
        new(a, A, z, B)
    end
end

struct TaylorInverter{N}
    f::Function
end

function TaylorInverter{N}() where {N}
    f = create_expressions(N)
    return TaylorInverter{N}(f)
end

function invert(ti::TaylorInverter, a::Vector)
    return ti.f([a])
end

function invert(ti::TaylorInverter, taylor1::Taylor1)
    order = taylor1.order
    a = taylor1.coeffs[2:end]
    substitution = Taylor1([-taylor1.coeffs[begin], 1], order)
    return Taylor1([0; ti.f([a])], order)(substitution)
end

function truncaterule(n, z)
    r = @acrule (~~a + ~b * (~z)^(~n::(m -> m > n))) => sum(~~a)
    rf = @acrule (~b * (~z)^(~n::(m -> m > n))) => 0
    return Fixpoint(Prewalk(Chain([r, rf])))
end

function truncatedpower(it::InverseTaylor, k::Int, truncrule::Fixpoint)
    single = mapreduce(kAn -> kAn[2] * it.z^kAn[1], +, enumerate(it.A))
    result = single
    for _ in 2:k
        result = simplify(
            result * single |> expand;
            rewriter=truncrule
        )
    end
    return result
end

function process(k::Int, an::Num, it::InverseTaylor, truncrule::Fixpoint)
    t = an * truncatedpower(it, k, truncrule)
    return simplify(
        t |> expand;
        rewriter=truncrule
    )
end

function initial_substitution(it::InverseTaylor{N}) where {N}
    truncrule = truncaterule(N, it.z)
    subbed = mapreduce(kan -> process(kan..., it, truncrule), +, enumerate(it.a)) |> expand

    subz = substitute(subbed, Dict(it.z => 0))
    subz == 0 || error("Function should be at least linear in z")

    return subbed
end

function further_substitution(it::InverseTaylor{N}, subbed::Num) where {N}
    for i in 1:N
        @info "Expanding term $i to create ivnersion expression"
        divided = simplify(subbed / it.z)
        subz = substitute(divided, Dict(it.z => 0))
        it.B[i] = solve_for(subz ~ i == 1 ? 1 : 0, it.A[i]) |> simplify

        subbed = substitute(divided - subz, Dict(it.A[i] => it.B[i]))
    end
    return subbed
end

function create_expressions(n::Int)
    it = InverseTaylor{n}()
    subbed = initial_substitution(it)
    subbed = further_substitution(it, subbed)

    f = build_function(it.B, [it.a])[1]  # [1] because build function also provided an in-place version in [2]
    return eval(f)
end

export TaylorInverter, invert

end # module TaylorInversion
