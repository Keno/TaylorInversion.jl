using ReTest

using TaylorInversion

include("TestTaylorInversion.jl")

retest(TaylorInversion, TestTaylorInversion)