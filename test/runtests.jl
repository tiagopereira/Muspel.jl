using Muspel
using Test

@testset "Muspel" verbose=true begin
    #include("read_utils.jl")
    #include("types.jl")
    #include("lte.jl")
    include("background.jl")
end
