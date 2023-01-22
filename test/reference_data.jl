function reference(a_trunc::Vector{<:Number})
    n = length(a_trunc)
    A = zeros(eltype(a_trunc), 7)
    a = zeros(eltype(a_trunc), 7)
    a[1:n] .= a_trunc
    A[1] = a[1]^(-1)
    A[2] = -a[1]^(-3) * a[2]
    A[3] = a[1]^(-5) * (2a[2]^2 - a[1]a[3])
    A[4] = a[1]^(-7) * (5a[1]a[2]a[3] - a[1]^2 * a[4] - 5a[2]^3)
    A[5] = a[1]^(-9) * (6a[1]^2 * a[2]a[4] + 3a[1]^2 * a[3]^2 + 14a[2]^4 - a[1]^3 * a[5] - 21a[1]a[2]^2 * a[3])
    A[6] = a[1]^(-11) * (7a[1]^3 * a[2]a[5] + 7a[1]^3 * a[3]a[4] + 84a[1]a[2]^3 * a[3] - a[1]^4 * a[6] - 28a[1]^2 * a[2]a[3]^2 - 42a[2]^5 - 28a[1]^2 * a[2]^2 * a[4])
    A[7] = a[1]^(-13) * (8a[1]^4 * a[2]a[6] + 8a[1]^4 * a[3]a[5] + 4a[1]^4 * a[4]^2 + 120a[1]^2 * a[2]^3 * a[4] + 180a[1]^2 * a[2]^2 * a[3]^2 + 132a[2]^6 - a[1]^5 * a[7] - 36a[1]^3 * a[2]^2 * a[5] - 72a[1]^3 * a[2]a[3]a[4] - 12a[1]^3 * a[3]^3 - 330a[1]a[2]^4 * a[3])
    return A[1:n]
end

function reference(a::Taylor1)
    A = reference(a.coeffs[2:end])
    order = length(A)
    return Taylor1([0; A], order)
end

taylor_inverters = Dict{Int64,TaylorInverter}()
function get_taylor_inverter(n::Int)
    return get(taylor_inverters, n, TaylorInverter{n}())
end