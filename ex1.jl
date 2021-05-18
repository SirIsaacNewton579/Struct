using LinearAlgebra
using PyPlot
PyPlot.pygui(true)
include("Structdefine.jl")
include("Structsolve.jl")

#ex6-5(a)
node = [[0.0;0],[2;0],[4,0],[1,1],[3,1]]
ecodes = [[1,2],[2,3],[1,4],[4,5],[2,4],[2,5],[3,5]]
el = build_truss(ecodes)
sp = [Support(1,[true,true,false]),Support(2,[false,true,false]),Support(3,[false,true,false])]
load=Load[Nodeforce(4,[0;-0.1;0]),Nodeforce(5,[0;-0.1;0])]
st = Struct(node,el,load,sp)

spforce,endpointforce,endpointdisplacement = solvestruct(st)
PyPlot.figure();plot_orign(st);plot_disp(st)
PyPlot.figure();plot_orign(st);plot_axis(st)
#ex6-5(b)
node = [[0.0;0],[1;0],[1;1],[0;1]]
ecodes = [[1,2],[2,3],[4,3],[4,1],[1,3],[2,4]]
el = build_truss(ecodes)
sp = [Support(1,[true,true,false]),Support(2,[false,true,false])]
load=Load[Nodeforce(4,[-0.1;0;0]),Nodeforce(3,[0.1;0;0])]
st = Struct(node,el,load,sp)
spforce,endpointforce,endpointdisplacement = solvestruct(st)
N,u,Q,M,θ,w = internalforcedispsingle(st,4)
PyPlot.figure();plot_orign(st);plot_disp(st)
PyPlot.figure();plot_orign(st);plot_axis(st)
#ex6-1(a)
node = [[0.0,0],[1,0]]
el = build_frame([[1,2]])
sp = [Support(1,[true,true,true]),Support(2,[false,true,false])]
load=[Concentratedforce(1,0.5,[0;-1;0],isgb)]
st = Struct(node,el,load,sp)
spforce,endpointforce,endpointdisplacement = solvestruct(st)
N,u,Q,M,θ,w = internalforcedispsingle(st,1)
PyPlot.figure();plot_orign(st);plot_shear(st)
PyPlot.figure();plot_orign(st);plot_moment(st)
#ex6-1(b)
k=2;a=1/3
node = [[0.0,0],[1-a,0],[1,0]]
el = build_frame([[1,2],[2,3]],[1.0,1],[1e12,1e12],[k,1.0])
sp = [Support(1,[true,true,true]),Support(3,[false,true,false])]
load=[Nodeforce(2,[0;-1;0])]
st = Struct(node,el,load,sp)
spforce,endpointforce,endpointdisplacement = solvestruct(st)
N,u,Q,M,θ,w = internalforcedispsingle(st,1)
PyPlot.figure();plot_orign(st);plot_shear(st)
PyPlot.figure();plot_orign(st);plot_moment(st)
#ex6-1(c)
node = [[0.0,0],[1,0]]
el = build_frame([[1,2]])
sp = [Support(1,[true,true,true]),Support(2,[true,true,true])]
load=[Concentratedforce(1,0.5,[0;-1;0],isgb)]
st = Struct(node,el,load,sp)
spforce,endpointforce,endpointdisplacement = solvestruct(st)
N,u,Q,M,θ,w = internalforcedispsingle(st,1)
PyPlot.figure();plot_orign(st);plot_shear(st)
PyPlot.figure();plot_orign(st);plot_moment(st)
#ex6-1(d)
node = [[0.0,0],[1,0]]
el = build_frame([[1,2]])
sp = [Support(1,[true,true,true]),Support(2,[true,false,true])]
load=[Concentratedforce(1,1,[0;-1;0],isgb)]
st = Struct(node,el,load,sp)
spforce,endpointforce,endpointdisplacement = solvestruct(st)
N,u,Q,M,θ,w = internalforcedispsingle(st,1)
PyPlot.figure();plot_orign(st);plot_shear(st)
PyPlot.figure();plot_orign(st);plot_moment(st)
#ex6-3(a)
node = [[0.0,0],[12,0],[0,6],[12,6]]
enodes = [[1,3],[3,4],[2,4]]
el = build_frame(enodes)
sp = [Support(1,[true,true,false]),Support(2,[true,true,false])]
load = [Uniformforce(1,[0,1],[0.1,0,0],isgb),Uniformforce(3,[0,1],[0.1,0,0],isgb)]
st = Struct(node,el,load,sp)
PyPlot.figure();plot_orign(st);plot_moment(st)
#ex6-3(b)
node = [[0.0,0],[1,0],[2,0],[1,-1]]
E = fill(1.0,3);A=fill(10.0,3);I=fill(1.0,3)
el = [build_frame([[1,2],[2,3]],E,A,I);build_truss([[2,4]],E,A,I)]
sp = [Support(1,[true,true,false]),Support(4,[true,true,false]),Support(3,[false,true,false])]
load = [Nodeforce(2,[0,-1,0])]
st = Struct(node,el,load,sp)
spforce,endpointforce,endpointdisplacement = solvestruct(st)
PyPlot.figure();plot_orign(st);plot_moment(st)
#ex6-3(c)
node = [[0.0,0],[15,0],[0,10],[2.5,10],[12.5,10],[15,10]]
enodes = [[1,4],[2,5],[3,4],[4,5],[5,6]]
el = build_frame(enodes)
sp = [Support(1,[true,true,false]),Support(2,[true,true,false])]
load = [Uniformforce(3,[0,1],[0,-0.01,0],isgb),Uniformforce(4,[0,1],[0,-0.01,0],isgb),Uniformforce(5,[0,1],[0,-0.01,0],isgb)]
st = Struct(node,el,load,sp)
spforce,endpointforce,endpointdisplacement = solvestruct(st)
N,u,Q,M,θ,w = internalforcedispsingle(st,1)
PyPlot.figure();plot_orign(st);plot_moment(st)
#ex6-4(a)
node = [[0.0,0],[12,0],[0,6],[12,6]]
el = [build_frame([[1,3],[2,4]]);build_truss([[3,4]],[1.0],[1e12],[1.0])]
sp = [Support(1,[true,true,true]),Support(2,[true,true,true])]
load = [Uniformforce(1,[0,1],[0.2,0,0],isgb)]
st = Struct(node,el,load,sp)
PyPlot.figure();plot_orign(st);plot_moment(st)
#ex6-4(b)
node = [[0.0,0.0],[12,0],[0,6],[12,6],[0,9],[12,9]]
enodes = [[1,3],[2,4],[3,5],[4,6],[5,6]];E=fill(1.0,5);A=fill(1e12,5);I=[2.0;2.0;fill(1.0,3)]
el = [build_frame(enodes[1:4],E[1:4],A[1:4],I[1:4]);build_truss([enodes[end]],[E[end]],[A[end]],[I[end]])]
sp = [Support(1,[true,true,true]),Support(2,[true,true,true])]
load = [Concentratedforce(4,1/3,[0.1,0,0],isgb)]
st = Struct(node,el,load,sp)
N,u,Q,M,θ,w=internalforcedispsingle(st,4)
PyPlot.figure();plot_orign(st);plot_moment(st)
#ex6-4(c)
node = [[0,0],[4,0],[10,0],[0,4.65],[0,6.75],[4,6.75],[10,6.75],[4,9.35],[10,9.35]]
enodes = [[1,4],[4,5],[5,6],[2,6],[6,8],[8,9],[3,7],[7,9]]
E = fill(1.0,8);A=fill(Inf,8);I=[2.83,1,1,8.10,1.59,1,8.10,1.59]
el = build_frame(enodes,E,A,I)
tsfidx=[3,6]
totruss!(el,tsfidx)
sp = [Support(1,[true,true,true]),Support(2,[true,true,true]),Support(3,[true,true,true])]
load = [Uniformforce(7,[0,1],[-1,0,0],isgb),Uniformforce(8,[0,1],[-1,0,0],isgb)]
st = Struct(node,el,load,sp)
spforce,endpointforce,endpointdisplacement = solvestruct(st)
N,u,Q,M,θ,w=internalforcedispsingle(st,7)
PyPlot.figure();plot_orign(st);plot_moment(st)