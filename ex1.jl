using LinearAlgebra
using PyPlot

includet("Structdefine.jl")
includet("Structsolve.jl")

#ex6-5(b)
node = [[0.0;0],[1;0],[1;1],[0;1]]
ecodes = [[1,2],[2,3],[4,3],[4,1],[1,3],[2,4]]
el = build_truss(ecodes)
sp = [Support(1,[true,true,false]),Support(2,[false,true,false])]
load=Load[Nodeforce(4,[-0.1;0;0]),Nodeforce(3,[0.1;0;0])]
st = Struct(node,el,load,sp)

spforce,endpointforce,endpointdisplacement = solvestruct(st)
N,u,Q,M,θ,w = internalforcedispsingle(st,4)
PyPlot.pygui(true)
plot_orign(st)
plot_disp(st)
PyPlot.figure()
plot_orign(st)
plot_axis(st)

#ex6-1(d)
node = [[0.0,0],[1,0]]
el = build_frame([[1,2]])
sp = [Support(1,[true,true,true]),Support(2,[true,false,true])]
load=[Concentratedforce(1,1,[0;-1;0],isgb)]
st = Struct(node,el,load,sp)
spforce,endpointforce,endpointdisplacement = solvestruct(st)
N,u,Q,M,θ,w = internalforcedispsingle(st,1)
PyPlot.figure()
plot_orign(st)
plot_shear(st)
PyPlot.figure()
plot_orign(st)
plot_moment(st)
#ex6-3(a)
node = [[0.0,0],[12,0],[0,6],[12,6]]
enodes = [[1,3],[3,4],[2,4]]
el = build_frame(enodes)
sp = [Support(1,[true,true,false]),Support(2,[true,true,false])]
load = [Uniformforce(1,[0,1],[0.1,0,0],isgb),Uniformforce(3,[0,1],[0.1,0,0],isgb)]
st = Struct(node,el,load,sp)
PyPlot.figure()
plot_orign(st)
plot_moment(st)
#ex6-3(b)
node = [[0.0,0],[1,0],[2,0],[1,-1]]
E = fill(1.0,3);A=fill(10.0,3);I=fill(1.0,3)
el = [build_frame([[1,2],[2,3]],E,A,I);build_truss([[2,4]],E,A,I)]
sp = [Support(1,[true,true,false]),Support(4,[true,true,false]),Support(3,[false,true,false])]
load = [Nodeforce(2,[0,-1,0])]
st = Struct(node,el,load,sp)
spforce,endpointforce,endpointdisplacement = solvestruct(st)
PyPlot.figure()
plot_orign(st)
plot_moment(st)
#ex6-3(c)
node = [[0.0,0],[15,0],[0,10],[2.5,10],[12.5,10],[15,10]]
enodes = [[1,4],[2,5],[3,4],[4,5],[5,6]]
el = build_frame(enodes)
sp = [Support(1,[true,true,false]),Support(2,[true,true,false])]
load = [Uniformforce(3,[0,1],[0,-0.01,0],isgb),Uniformforce(4,[0,1],[0,-0.01,0],isgb),Uniformforce(5,[0,1],[0,-0.01,0],isgb)]
st = Struct(node,el,load,sp)
spforce,endpointforce,endpointdisplacement = solvestruct(st)
N,u,Q,M,θ,w = internalforcedispsingle(st,1)
PyPlot.figure()
plot_orign(st)
plot_moment(st)