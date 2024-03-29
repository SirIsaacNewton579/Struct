
struct Element
    enode :: Vector{Int64}
    E :: Float64
    A :: Float64
    I :: Float64
    constraint :: Vector{Vector{Bool}}
end

isgb = true;islc = false
abstract type Load end
struct Nodeforce <: Load
    fnode :: Int64
    magnitude :: Vector{Float64}
end
abstract type ElementLoad <: Load end
struct Concentratedforce <: ElementLoad
    elementcode :: Int64
    loaction :: Float64
    magnitude :: Vector{Float64}
    fgol :: Bool
end
function Concentratedforce(elc,loac,mag;fgol=isgb)
    Concentratedforce(elc,loac,mag,fgol)
end

struct Disturbutionforce <: ElementLoad
    elementcode :: Int64
    loaction :: Vector{Float64}
    magnitude :: Function
    fgol :: Bool
end
function Disturbutionforce(elc,loac,mag;fgol=isgb)
    Disturbutionforce(elc,loac,mag,fgol)
end
struct Uniformforce <: ElementLoad
    elementcode :: Int64
    loaction :: Vector{Float64}
    magnitude :: Vector{Float64}
    fgol :: Bool
end
function Uniformforce(elc,loac,mag;fgol=isgb)
    Uniformforce(elc,loac,mag,fgol)
end
struct Support
    snode :: Int64
    ncst :: Vector{Bool}
    angle :: Float64
end
function Support(sn,nc;ag=0.0)
    Support(sn,nc,ag)
end

struct Struct
    globalcoordinate :: Vector{Vector{Float64}}
    element :: Vector{Element}
    load :: Vector{Load}
    support :: Vector{Support}
end
function build_truss(enodes::Vector{Vector{Int64}})
    elnum = length(enodes)
    E = fill(1.0,elnum)
    A = fill(1.0,elnum)
    I = fill(1.0,elnum)
    build_truss(enodes,E,A,I)
end
function build_truss(enodes::Vector{Vector{Int64}},E::Vector{Float64},A::Vector{Float64},I::Vector{Float64})
    elnum = length(enodes)
    el = Vector{Element}(undef,elnum)
    constr = [[true,true,false],[true,true,false]]
    for i = 1:elnum
        el[i]=Element(enodes[i],E[i],A[i],I[i],constr)
    end
    el
end
function build_frame(enodes::Vector{Vector{Int64}})
    elnum = length(enodes)
    E = fill(1.0,elnum)
    A = fill(1e12,elnum)
    I = fill(1.0,elnum)
    build_frame(enodes,E,A,I)
end
function build_frame(enodes::Vector{Vector{Int64}},E::Vector{Float64},A::Vector{Float64},I::Vector{Float64})
    elnum = length(enodes)
    el = Vector{Element}(undef,elnum)
    constr = [[true,true,true],[true,true,true]]
    for i = 1:elnum
        el[i]=Element(enodes[i],E[i],A[i],I[i],constr)
    end
    el
end
function totruss!(el,tsfidx)
    for i in tsfidx
        el[i] = Element(el[i].enode,el[i].E,el[i].A,el[i].I,[[true,true,false],[true,true,false]])
    end
end
function build_fix(sn::Vector{Int64})
    sp=Vector{Support}(undef,length(sn))
    for i=1:length(sn)
        sp[i]=Support(sn[i],fill(true,3))
    end
    sp
end
function build_hinge(sn::Vector{Int64})
    sp=Vector{Support}(undef,length(sn))
    for i=1:length(sn)
        sp[i]=Support(sn[i],[true,true,false])
    end
    sp
end