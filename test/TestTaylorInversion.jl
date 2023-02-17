module TestTaylorInversion
using ReTest
using TaylorSeries
import Base.isapprox

using TaylorInversion

include("reference_data.jl")

include("test_inversion.jl")
include("test_taylor1.jl")

end  # module TestTaylorInversion
