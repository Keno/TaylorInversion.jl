module TestTaylorInversion

using TaylorSeries

function reference(a::Vector{<:Number})
    A[1] = a[1]
    A[2] = -a[1]^(-3)a[2] 
    A[3] = a[1]^(-5)(2a[2]^2-a[1]a[3]) 
    A[4] = a[1]^(-7)(5a[1]a[2]a[3]-a[1]^2a[4]-5a[2]^3) 
    A[5] = a[1]^(-9)(6a[1]^2a[2]a[4]+3a[1]^2a[3]^2+14a[2]^4-a[1]^3a[5]-21a[1]a[2]^2a[3]) 
    A[6] = a[1]^(-11)(7a[1]^3a[2]a[5]+7a[1]^3a[3]a[4]+84a[1]a[2]^3a[3]-a[1]^4a[6]-28a[1]^2a[2]a[3]^2-42a[2]^5-28a[1]^2a[2]^2a[4]) 
    A[7] = a[1]^(-13)(8a[1]^4a[2]a[6]+8a[1]^4a[3]a[5]+4a[1]^4a[4]^2+120a[1]^2a[2]^3a[4]+180a[1]^2a[2]^2a[3]^2+132a[2]^6-a[1]^5a[7]-36a[1]^3a[2]^2a[5]-72a[1]^3a[2]a[3]a[4]-12a[1]^3a[3]^3-330a[1]a[2]^4a[3]) 
    return A
end

function reference(a::Taylor1)
    A = reference(a.coeffs[2:end])
    order = length(A)
    return Taylor1([0; A], order)
end

n = 5
@variables x, y, a[1:n], A[1:n]
B = zeros(Num, n)

y = mapreduce(kan -> kan[2] * x^kan[1], +, enumerate(a))
subbed = substitute(y, Dict(x=>mapreduce(kAn -> kAn[2] * z^kAn[1], +, enumerate(A))))


this = nothing
remaining2 = nothing
subbed2 = nothing
divided = nothing
subtracted = nothing

rule1 = @acrule (~x + ~y)/(~y) => (~x)/(~y) + 1
rule2 = @acrule (~x + ~y + ~z)/(~y) => (~x + ~z)/(~y) + 1
cr = SymbolicUtils.Chain([rule1, rule2])

rule3 = @acrule (~x + ~y)/(~z) => ~x / ~z + ~y / ~z

remaining = subbed
for i in 1:n
    differential = Symbolics.Differential(z)(remaining) |> expand_derivatives
    this = substitute(differential, z=>0)
    target = i == 1 ? 1 : 0
    B[i] = Symbolics.solve_for(this ~ target, A[i]) |> simplify
    subbed2 = substitute(remaining, Dict(A[i] => B[i])) |> simplify
    subtracted = expand(subbed2) - z
    divided = expand(subbed2) / z |> simplify
    remaining = divided
end

end  # module TestTaylorInversion
