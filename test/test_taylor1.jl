function rand(::Type{Taylor1}, n::Int)
    a = randn(n + 1)
    return Taylor1(a, n)
end

function isapprox(left::Taylor1, right::Taylor1; kwargs...)
    return left.order == right.order && isapprox(left.coeffs, right.coeffs; kwargs...)
end

n = 4
@time ti = TaylorInversion.TaylorInverter{n}()
@testset verbose = false "Compare Taylor1 to reference" begin
    @testset "Random number $i_exp" for i_exp in 1:10
        taylor1 = rand(Taylor1, n)
        @test isapprox(invert(ti, taylor1), reference(taylor1))
    end
end

@testset verbose = false "Test Taylor1 self-consistency" begin
    @testset "Random number $i_exp" for i_exp in 1:1
        taylor1 = rand(Taylor1, n)
        @test isapprox(invert(ti, invert(ti, taylor1)), taylor1)
    end
end