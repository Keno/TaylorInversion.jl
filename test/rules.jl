using SymbolicUtils
using SymbolicUtils: Chain, RestartedChain, Fixpoint, Prewalk

@syms x::Real y::Real z::Real


@syms a b c d e f g

function getrule(m::Int)
    r1 = @acrule ~b * x^m + ~c * x^m => (~b + ~c) * x^m
    # r1 = @acrule ~~a1 + ~b * x^m + ~~a2 + ~c * x^m => sum(~~a1) + (~b + ~c) * x^m;
    # r2 = @acrule ~~a1 + ~b * x^m + ~~a2 + ~c * x^m + ~~a3::(y->length(y) > 0) => sum(~~a1) + (~b + ~c) * x^m + sum(~~a3);
    r3 = @acrule ~~a + ~b * x^m +  ~  ~ c => ~b * x^x
    return Fixpoint(RestartedChain([r1])), r3
end

# fp, r = getrule(n)
rfp = @acrule ~b * x^~n + ~c * x^~n => (~b + ~c) * x^~n;
rfp1 = @acrule ~b * x + ~c * x => (~b + ~c) * x;
fp = Fixpoint(Chain([rfp, rfp1]))

@show fp(a + b * x^2 + c * x^1 + d)
@show fp(b * x^1 + c + d * x^1)
@show fp(b * x^2 + c + d * x^2 + e + f)
@show fp(b * x^2 + c + b)
@show fp(a + b * x^3 + c * x^2 + d + e + f * x^3 + g * x^2)

rnb = @acrule (~~a + ~b * x^~n) / x => sum(~~a) / x + ~b * x^(~n - 1)
rn1 = @acrule (~~a + x^~n) / x => sum(~~a) / x + x^(~n - 1)
r1b = @acrule (~~a + ~b * x) / x => sum(~~a) / x + ~b
r11 = @acrule (~~a + x) / x => sum(~~a) / x + 1
fp = Fixpoint(Prewalk(Chain([rnb, rn1, r1b, r11])))
@show fp((a * x^2 + x + b * x^3) / x + 1 + (b * x^2 + b * x) / x)