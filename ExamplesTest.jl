

include("2D\ Cauchy\ Distribution.jl")
include("Airy\ equation.jl")
if isdir(Pkg.dir("RandomMatrices"))
    include("GUE\ Sampling.jl")
end
include("Lanczos.jl")
include("Lee\ &\ Greengard\ equation.jl")
include("Newton iteration for Nonlinear BVPs.jl")
