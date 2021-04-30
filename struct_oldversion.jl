


struct Element
    enode :: Vector{Int64}
    E :: Float64
    A :: Float64
    I :: Float64
    constraint :: Vector{Vector{Bool}}
end

abstract type Load end

struct Concentratedforce <: Load
    elementcode :: Int64
    loaction :: Float64
    magnitude :: Vector{Float64}
end


struct Node
    state :: Vector{Bool}
    coordinate :: Vector{Float64}
end

struct Struct
    globalnode :: Vector{Node}
    element :: Vector{Element}
    load :: Vector{Load}
end

