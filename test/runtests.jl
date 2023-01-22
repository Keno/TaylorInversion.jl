using ReTest

using TaylorInversion
# (@__DIR__) in LOAD_PATH || push!(LOAD_PATH, @__DIR__)
# @show LOAD_PATH
include("TestTaylorInversion.jl")
# using TestTaylorInversion

retest(TaylorInversion, TestTaylorInversion)