@testset "Test object creation for n=$n" for n in 1:3
    @test isa(InverseTaylor{n}(), InverseTaylor)
end

@testset "Test for n=$n" for n in 1:7
    @time f = TaylorInversion.create_expressions(n)
    @testset verbose = false "Testing various random values for a" begin
        @testset "Random number $i_exp" for i_exp in 1:10
            a = randn(n)
            @test isapprox(f([a]), reference(a))
        end
    end
end