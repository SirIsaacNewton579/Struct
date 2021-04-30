
struct Element
    enode :: Vector{Int64}
    E :: Float64
    A :: Float64
    I :: Float64
    constraint :: Vector{Vector{Bool}}
end

isgb = true;islc = false
abstract type Load end

struct Concentratedforce <: Load
    elementcode :: Int64
    loaction :: Float64
    magnitude :: Vector{Float64}
    fgol :: Bool
end


struct Disturbutionforce <: Load
    elementcode :: Int64
    loaction :: Vector{Float64}
    magnitude :: Function
    fgol :: Bool
end
struct Uniformforce <: Load
    elementcode :: Int64
    loaction :: Vector{Float64}
    magnitude :: Vector{Float64}
    fgol :: Bool
end

struct Support
    snode :: Int64
    ncst :: Vector{Bool}
end

struct Struct
    globalcoordinate :: Vector{Vector{Float64}}
    element :: Vector{Element}
    load :: Vector{Load}
    support :: Vector{Support}
end