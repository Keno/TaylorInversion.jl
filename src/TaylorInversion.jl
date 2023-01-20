module TaylorInversion
using Symbolics: @variables, @acrule, @rule
using Symbolics: Num, Chain, Fixpoint, Prewalk
using Symbolics: build_function, solve_for, expand, simplify, substitute

function truncaterule(n, z)
    r = @acrule (~~a + ~b * (~z)^(~n::(m -> m > n))) => sum(~~a)
    rf = @acrule (~b * (~z)^(~n::(m -> m > n))) => 0
    return Fixpoint(Prewalk(Chain([r, rf])))
end

function process(kan, A, n, z)
    k, an = kan
    t = an * mapreduce(kAn -> kAn[2] * z^kAn[1], +, enumerate(A))^k
    return simplify(
        t |> expand;
        rewriter=truncaterule(n, z)
    )
end

function create_expressions(n::Int)
    @variables x, y, z, a[1:n], A[1:n]
    B = zeros(Num, n)

    subbed = mapreduce(kan -> process(kan, A, n, z), +, enumerate(a)) |> expand

    subz = substitute(subbed, Dict(z => 0))
    subz == 0 || error("Function should be at least linear in z")

    # fp = getfp()
    for i in 1:n
        @info "Expanding term $i to create ivnersion expression"
        divided = simplify(subbed / z)
        subz = substitute(divided, Dict(z => 0))
        B[i] = solve_for(subz ~ i == 1 ? 1 : 0, A[i]) |> simplify

        subbed = substitute(divided - subz, Dict(A[i] => B[i]))
    end

    f = build_function(B, [a])[1] |> eval
    @show B
    return f
end

end # module TaylorInversion
