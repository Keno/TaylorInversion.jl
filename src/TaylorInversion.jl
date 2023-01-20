module TaylorInversion
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

function truncaterule(n, z)
    r = @acrule (~~a + ~b * (~z)^(~n::(m -> m > n))) => sum(~~a)
    rf = @acrule (~b * (~z)^(~n::(m -> m > n))) => 0
    return Fixpoint(Prewalk(Chain([r, rf])))
end

function truncatedpower(A, z, k::Int, truncrule)
    single = mapreduce(kAn -> kAn[2] * z^kAn[1], +, enumerate(A))
    result = single
    for kk in 2:k
        result = simplify(
            result * single |> expand;
            rewriter=truncrule
        )
    end
    return result
end

function process(kan, A, truncrule, z)
    k, an = kan
    t = an * truncatedpower(A, z, k, truncrule)
    return simplify(
        t |> expand;
        rewriter=truncrule
    )
end

function initial_substitution(n::Int, a::Arr{Num,1}, A::Arr{Num,1}, z::Num)
    truncrule = truncaterule(n, z)
    subbed = mapreduce(kan -> process(kan, A, truncrule, z), +, enumerate(a)) |> expand

    subz = substitute(subbed, Dict(z => 0))
    subz == 0 || error("Function should be at least linear in z")

    return subbed
end

function create_expressions(n::Int)
    @variables z, a[1:n], A[1:n]
    subbed = initial_substitution(n, a, A, z)

    B = zeros(Num, n)
    for i in 1:n
        @info "Expanding term $i to create ivnersion expression"
        divided = simplify(subbed / z)
        subz = substitute(divided, Dict(z => 0))
        B[i] = solve_for(subz ~ i == 1 ? 1 : 0, A[i]) |> simplify

        subbed = substitute(divided - subz, Dict(A[i] => B[i]))
    end

    f = build_function(B, [a])[1] |> eval
    return f
end

export InverseTaylor

end # module TaylorInversion
