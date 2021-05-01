using LinearAlgebra
using PyPlot

include("Structdefine.jl")
include("Structsolve.jl")


st = Struct(
    [[0.0,0.0],[1.0,0.0],[0.0,1.0],[1.0,1.0]],
    [Element([1,3],1e9,1e-3,1e-6,[[true,true,false],[true,true,true]]),
    Element([2,4],1e9,1e-3,1e-6,[[true,true,false],[true,true,true]]),
    Element([3,4],1e9,1e-3,1e-6,[[true,true,true],[true,true,true]]),
    Element([1,4],1e9,1e-3,1e-6,[[true,true,false],[true,true,false]]),
    Element([2,3],1e9,1e-3,1e-6,[[true,true,false],[true,true,false]])],
    [Disturbutionforce(3,[0,1],s->[100e3+0*s,0,0],isgb)],
    [Support(1,[true,true,false]),Support(2,[true,true,false])]
)

node = st.globalcoordinate;el = st.element;load = st.load
sp = st.support
spb = [sup.ncst for sup in sp]
spbindex = findall.(spb)
svalue = length.(spbindex)
nvalue = 12*length(el) + sum(svalue)

function build_A(st :: Struct)
A = Vector{Float64}[]
#位置量顺序：支座力，[u1,v1,θ1,fx1,fy1,M1,u2,v2,θ2,fx2,fy2,M2]
node = st.globalcoordinate;el = st.element;load = st.load
sp = st.support
spb = [sup.ncst for sup in sp]
spbindex = findall.(spb)
svalue = length.(spbindex)
nvalue = 12*length(el) + sum(svalue)
build_A_consider_support!(A,sp,el)
build_A_consider_connect!(A,node,el,sp,nvalue,svalue)
build_A_consider_equilibrium_and_displacement!(A,node,el,load,nvalue,svalue)
A
end

A = build_A(st)
Ac = transpose(reduce(hcat,A))
B = Ac[:,1:end-1]
b = Ac[:,end]

x = B\b
ssval = sum(svalue)
nodeforce = x[1:ssval]
se = [x[ssval+1+12(i-1):ssval+12+12(i-1)] for i =1:length(el)]


node_index = [nodeonele(i,el) for i = 1:length(node)]
u = [se[findfirst(!isempty,ni)][(1:2) .+ 6*(first(ni[findfirst(!isempty,ni)])-1)] for ni in node_index]
PyPlot.pygui(true)
PyPlot.figure()
for eli in el
    plot([node[eli.enode[1]][1],node[eli.enode[2]][1]],[node[eli.enode[1]][2],node[eli.enode[2]][2]])
end
for eli in el
    nnode = node + u
    plot([nnode[eli.enode[1]][1],nnode[eli.enode[2]][1]],[nnode[eli.enode[1]][2],nnode[eli.enode[2]][2]])
end