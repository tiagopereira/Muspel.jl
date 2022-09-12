using Muspel
using Test

@testset "Muspel" verbose=true begin
    include("types.jl")
    include("read_utils.jl")
    include("lte.jl")
    include("background.jl")
end
