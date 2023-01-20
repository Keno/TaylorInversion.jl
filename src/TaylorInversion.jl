module TaylorInversion
using Symbolics: @variables
using Symbolics: Num, build_function
using Symbolics: solve_for, expand, simplify, substitute

function create_expressions(n::Int)
    @variables x, y, z, a[1:n], A[1:n]
    B = zeros(Num, n)

    y = mapreduce(kan -> kan[2] * x^kan[1], +, enumerate(a))
    subbed = substitute(y, Dict(x => mapreduce(kAn -> kAn[2] * z^kAn[1], +, enumerate(A)))) |> expand

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
