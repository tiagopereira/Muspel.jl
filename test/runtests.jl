using Muspel
using Test

@testset "Muspel" verbose=true begin
    # Ordered by dependency -> Should fix top tests first
    #include("types.jl")
    include("read_utils.jl")
    include("lte.jl")
    include("background.jl")
end
