@testset "Truncatedpower time for n=$n" for n in 8:9
    @testset "k = $k" for k in n-1:n
        it = TaylorInversion.InverseTaylor{n}()
        truncrule = TaylorInversion.truncaterule(n, it.z)
        @time segmented = TaylorInversion.truncatedpower(it, k, truncrule)
        @time baseline = TaylorInversion.truncatedpower(it, k, truncrule, Val{:baseline}())
        @test baseline - segmented == 0
    end
end

@testset "Test object creation for n=$n" for n in 1:3
    @test isa(TaylorInversion.InverseTaylor{n}(), TaylorInversion.InverseTaylor)
    # use the shared get_taylor_inverter to be able to reuse the objects later
    @time it = get_taylor_inverter(n)
    @test isa(it, TaylorInverter)
end


@testset "Test for n=$n" for n in 1:7
    ti = get_taylor_inverter(n)
    @testset verbose = false "Testing various random values for a" begin
        @testset "Random number $i_exp (n=$n)" for i_exp in 1:1
            a = randn(n)
            @test isapprox(invert(ti, a), reference(a))
        end
    end
end