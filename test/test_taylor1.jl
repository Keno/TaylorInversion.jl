function rand(::Type{Taylor1}, n::Int)
    a = randn(n + 1)
    a[1] = 0
    return Taylor1(a, n)
end

function Base.isapprox(left::Taylor1, right::Taylor1; kwargs...)
    return left.order == right.order && isapprox(left.coeffs, right.coeffs; kwargs...)
end

@testset verbose = false "Compare Taylor1 to reference" begin
    n = 3
    ti = get_taylor_inverter(n)
    @testset "Random number $i_exp" for i_exp in 1:10
        taylor1 = rand(Taylor1, n)
        @test isapprox(invert(ti, taylor1), reference(taylor1))
    end
end

@testset verbose = false "Test Taylor1 self-consistency" begin
    n = 3
    ti = get_taylor_inverter(n)
    @testset "Random number $i_exp" for i_exp in 1:10
        taylor1 = rand(Taylor1, n)
        taylor1.coeffs[1] = 0  # the inversion of the Taylor Series is agnostic of x0, y0.
        @test isapprox(invert(ti, invert(ti, taylor1)), taylor1)
    end
end

@testset verbose = false "Test Taylor1 for quadratics" begin
    n = 3
    ti = get_taylor_inverter(n)
    triplets = [(1, 1, 1), (1, 0, 1), (1, 0, 0), (1, 1, 0)]
    @testset "Test for [a, b, c] = [$a, $b, $c]" for (a, b, c) in triplets
        c, b, a = 0, 1, 1
        taylor1 = Taylor1([c, b, a, 0, 0])
        invtaylor1 = invert(ti, taylor1)
        refs = 0:0.01:0.1
        samples = taylor1.(refs)
        @test isapprox(invtaylor1.(samples), refs; rtol=1.e-2)
    end
end

@testset "Test Taylor1 linear" begin
    n = 3
    ti = get_taylor_inverter(n)
    c, b = 1, 1
    taylor1 = Taylor1([c, b], 5)
    invtaylor1 = invert(ti, taylor1)
    refs = 0:0.01:0.2
    samples = taylor1.(refs)
    @test isapprox(invtaylor1.(samples), refs; atol=1.e-15)
end