using LinearAlgebra
using PyPlot

include("Structdefine.jl")
include("Structsolve.jl")

node = [[0.0,0.0],[2.0,0.0],[0.0,1.0],[2.0,1.0]]
el = [Element([1,3],1e9,1e-3,1e-6,[[true,true,false],[true,true,true]]),
Element([2,4],1e9,1e-3,1e-6,[[true,true,false],[true,true,true]]),
Element([3,4],1e9,1e-3,1e-6,[[true,true,true],[true,true,true]]),
Element([1,4],1e9,1e-3,1e-6,[[true,true,false],[true,true,false]]),
Element([2,3],1e9,1e-3,1e-6,[[true,true,false],[true,true,false]])]
load = [Disturbutionforce(4,[0,1],s->[0,-0.5e3s,0],islc),Uniformforce(3,[0,1],[1e3,-1e3,0],islc)]
sp = [Support(1,[true,true,false]),Support(2,[true,true,false])]
st = Struct(node,el,load,sp)

spforce,endpointforce,endpointdisplacement = solvestruct(st)
N,u,Q,M,θ,w = internalforcedispsingle(node,el[4],4,load,endpointforce[4],endpointdisplacement[4])
PyPlot.pygui(true)
plot_orign(st)
plot_disp(st)
PyPlot.figure()
plot_orign(st)
plot_axis(st)
PyPlot.figure()
plot_orign(st)
plot_shear(st)
PyPlot.figure()
plot_orign(st)
plot_moment(st)

st2 = Struct([[0.0,0.0],[1.0,0.0]],
[Element([1,2],1e9,1e-3,1e-6,[[true,true,false],[false,true,false]])],
[Uniformforce(1,[0,1],[0,-1e3,0],isgb)],
[Support(1,[true,true,false]),Support(2,[false,true,false])]
)
build_A(st2)
solvestruct(st2)
plot_orign(st2)
plot_moment(st2)
plot_shear(st2)