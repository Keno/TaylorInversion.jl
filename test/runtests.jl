using ReTest

using TaylorInversion
@show push!(LOAD_PATH, @__DIR__)
using TestTaylorInversion

retest(TaylorInversion, TestTaylorInversion)