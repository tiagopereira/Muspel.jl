using RomanNumerals

# From DimensionalData.jl
function dimsymbols(i)
    symbols = ['↓', '→', '↗', '⬔', '◩', '⬒', '⬓', '■']
    symbols[min(i, length(symbols))]
end

function dimcolors(i)
    colors = [209, 32, 81, 204, 249, 166, 37]
    c = rem(i - 1, length(colors)) + 1
    colors[c]
end

function print_sizes(io, size;
    colors=map(dimcolors, ntuple(identity, length(size)))
)
    block_size = 1
    if !isempty(size)
        foreach(enumerate(size[1:end-1])) do (n, s)
            printstyled(io, s; color=colors[n])
            print(io, '×')
            block_size += length(string(s)) + 1
        end
        printstyled(io, last(size); color=colors[length(size)])
        block_size += length(string(last(size)))
    end
    return block_size
end

# From YAXArrays.jl
function formatbytes(x)
    exts = ["bytes", "KB", "MB", "GB", "TB","PB"]
    i = 1
    while x >= 1024
        i = i + 1
        x = x / 1024
    end
    return string(round(x, digits=2), " ", exts[i])
end

##
## Show functions for different types
##
function Base.show(io::IO, ::MIME"text/plain", a::A) where {A <: Muspel.AbstractAtmos}
    printstyled(io, "┌ "; color=:light_black)
    sum_width = print_sizes(io, (a.nz, a.ny, a.nx))
    summary_line = " $(string(nameof(typeof(a)))){$(eltype(a.temperature))}"
    print(io, summary_line)
    printstyled(io, " ┐"; color=:light_black)
    println(io)
    line_width = displaysize(io)[2]
    hdr_width = sum_width + length(summary_line) + 1
    top_line = string('├', '─'^max(0, hdr_width), '┴', '─'^max(0, line_width - hdr_width - 9), " dims ┐")
    printstyled(io, top_line; color=:light_black)
    println(io)

    # dimensions
    printstyled(io, "  ↓ z "; color=dimcolors(1))
    min_z = round(minimum(a.z) / 1e6, sigdigits=4)
    max_z = round(maximum(a.z) / 1e6, sigdigits=4)
    print(io, a.nz, " points: $min_z, …, $max_z Mm")
    println(io)

    printstyled(io, "  → y "; color=dimcolors(2))
    if a.ny == 1
        print(io, a.ny, " point")
        else
        if typeof(a) <: Atmosphere3D
            min_y = round(minimum(a.y) / 1e6, sigdigits=4)
            max_y = round(maximum(a.y) / 1e6, sigdigits=4)
            print(io, a.ny, " points: $min_y, …, $max_y Mm")
        else
            print(io, a.ny, " points")
        end
    end
    println(io)

    printstyled(io, "  ↗ x "; color=dimcolors(3))
    if a.nx == 1
        print(io, a.nx, " point")
    else
        if typeof(a) <: Atmosphere3D
            min_x = round(minimum(a.x) / 1e6, sigdigits=4)
            max_x = round(maximum(a.x) / 1e6, sigdigits=4)
            print(io, a.nx, " points: $min_x, …, $max_x Mm")
        else
            print(io, a.nx, " points")
        end
    end
    println(io)

    # memory usage
    mid_line = string('├', '─'^max(0, line_width -  20), " loaded in memory ┤")
    printstyled(io, mid_line; color=:light_black)
    println(io)
    print(io, "  data size: $(formatbytes(Base.summarysize(a)))")
    println(io)

    bottom_line = string('└', '─'^max(0, line_width - 2), '┘')
    printstyled(io, bottom_line; color=:light_black)
end


function Base.show(io::IO, ::MIME"text/plain", a::A) where {A <: AtomicModel}
    printstyled(io, "┌ "; color=:light_black)
    sum_width = print_sizes(io, (a.nlevels, a.nlines, a.ncontinua))
    summary_line = " $(string(nameof(typeof(a)))){$(eltype(a.χ)), $(eltype(a.g))}"
    print(io, summary_line)
    printstyled(io, " ┐"; color=:light_black)
    println(io)
    line_width = displaysize(io)[2]
    hdr_width = sum_width + length(summary_line) + 1
    top_line = string('├', '─'^max(0, hdr_width), '┴', '─'^max(0, line_width - hdr_width - 9), " dims ┐")
    printstyled(io, top_line; color=:light_black)
    println(io)

    stages = join(string.(RomanNumeral.(minimum(a.stage):maximum(a.stage))), ", ")
    print(io, "  Element: $(a.element) (", stages, ")")
    println(io)

    # dimensions
    printstyled(io, "  ⩸ ", a.nlevels, " levels:"; color=dimcolors(1))
    min_l = round(minimum(a.χ) * 1e18, sigdigits=3)
    max_l = round(maximum(a.χ) * 1e18, sigdigits=3)
    print(io, " $min_l, …, $max_l aJ")
    println(io)

    if !isempty(a.lines)
        printstyled(io, "  ↥ ", a.nlines, " lines:"; color=dimcolors(2))
        waves = [l.λ0 for l in a.lines]
        min_w = round(minimum(waves), sigdigits=4)
        max_w = round(maximum(waves), sigdigits=4)
        print(io, " $min_w, …, $max_w nm")
        println(io)
    end

    if !isempty(a.continua)
        printstyled(io, "  ⩘ ", a.ncontinua, " continua:"; color=dimcolors(3))
        waves = [l.λedge for l in a.continua]
        min_w = round(minimum(waves), sigdigits=4)
        max_w = round(maximum(waves), sigdigits=4)
        print(io, " $min_w, …, $max_w nm")
        println(io)
    end

    bottom_line = string('└', '─'^max(0, line_width - 2), '┘')
    printstyled(io, bottom_line; color=:light_black)
end


function Base.show(io::IO, ::MIME"text/plain", a::A) where {A <: AtomicLine}
    printstyled(io, "┌ "; color=:light_black)
    sum_width = 1
    summary_line = "$(string(nameof(typeof(a)))){$(eltype(a.λ))}"
    print(io, summary_line)
    printstyled(io, " ┐"; color=:light_black)
    println(io)
    line_width = displaysize(io)[2]
    hdr_width = sum_width + length(summary_line) + 1
    top_line = string('├', '─'^max(0, hdr_width), '┴', '─'^max(0, line_width - hdr_width - 13), " metadata ┐")
    printstyled(io, top_line; color=:light_black)
    println(io)

    # dimensions
    printstyled(io, "  ↥ ", "λ₀ at"; color=dimcolors(3))
    print(io, " $(round(a.λ0, sigdigits=4)) nm")
    nzeeman = max(1, length(a.σr_strength) + length(a.π_strength) + length(a.σb_strength))
    if nzeeman > 1
        print(io, ", $(nzeeman) Zeeman components")
    else
        print(io, ", $(nzeeman) Zeeman component")
    end
    println(io)

    printstyled(io, "  ∿ ", a.nλ, " wavelengths:"; color=dimcolors(1))
    min_w = round(minimum(a.λ), sigdigits=4)
    max_w = round(maximum(a.λ), sigdigits=4)
    print(io, " $min_w, …, $max_w nm")
    println(io)

    bottom_line = string('└', '─'^max(0, line_width - 2), '┘')
    printstyled(io, bottom_line; color=:light_black)
end


function Base.show(io::IO, ::MIME"text/plain", a::A) where {A <: AtomicContinuum}
    printstyled(io, "┌ "; color=:light_black)
    sum_width = 1
    summary_line = "$(string(nameof(typeof(a)))){$(eltype(a.σ))}"
    print(io, summary_line)
    printstyled(io, " ┐"; color=:light_black)
    println(io)
    line_width = displaysize(io)[2]
    hdr_width = sum_width + length(summary_line) + 1
    top_line = string('├', '─'^max(0, hdr_width), '┴', '─'^max(0, line_width - hdr_width - 13), " metadata ┐")
    printstyled(io, top_line; color=:light_black)
    println(io)

    # dimensions
    printstyled(io, "  ⩘ ", "edge at"; color=dimcolors(3))
    print(io, " $(round(a.λedge, sigdigits=4)) nm")
    println(io)

    printstyled(io, "  ∿ ", a.nλ, " wavelengths:"; color=dimcolors(1))
    min_w = round(minimum(a.λ), sigdigits=4)
    max_w = round(maximum(a.λ), sigdigits=4)
    print(io, " $min_w ... $max_w nm")
    println(io)

    bottom_line = string('└', '─'^max(0, line_width - 2), '┘')
    printstyled(io, bottom_line; color=:light_black)
end


function Base.show(
    io::IO,
    ::MIME"text/plain",
    v::AbstractVector{A},
) where {A <: Union{AtomicLine, AtomicContinuum}}
    nlines = length(v)
    nwave = sum([a.nλ for a in v])
    summary_line = "$nlines-element Vector{$(string(nameof(A)))} ($nwave wavelengths)"
    print(io, summary_line)
    println(io)
    if A <: AtomicLine
        println(io, " λ₀ (nm):")
        λ = round.([a.λ0 for a in v], sigdigits=6)
    elseif A <: AtomicContinuum
        println(io, " ⩘ λ (nm):")
        λ = round.([a.λedge for a in v], sigdigits=6)
    end
    Base.print_matrix(io, λ')
end
