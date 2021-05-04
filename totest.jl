using LinearAlgebra
using PyPlot

includet("Structdefine.jl")
includet("Structsolve.jl")


st = Struct(
    [[0.0,0.0],[2.0,0.0],[0.0,1.0],[2.0,1.0]],
    [Element([1,3],1e9,1e-3,1e-6,[[true,true,false],[true,true,true]]),
    Element([2,4],1e9,1e-3,1e-6,[[true,true,false],[true,true,true]]),
    Element([3,4],1e9,1e-3,1e-6,[[true,true,true],[true,true,true]]),
    Element([1,4],1e9,1e-3,1e-6,[[true,true,false],[true,true,false]]),
    Element([2,3],1e9,1e-3,1e-6,[[true,true,false],[true,true,false]])],
    [Disturbutionforce(4,[0,1],s->[0,-0.5e3s^2,0],islc),Disturbutionforce(3,[0,1],s->[0,-1e3+0*s,0],islc),Disturbutionforce(3,[0,1],s->[1e3+0*s,0,0],islc)],
    [Support(1,[true,true,false]),Support(2,[true,true,false])]
)

node = st.globalcoordinate;el = st.element;load = st.load
sp = st.support
spforce,endpointforce,endpointdisplacement = solvestruct(st)
N,u,Q,M,θ,w = internalforcedispsingle(node,el[3],3,load,endpointforce[3],endpointdisplacement[3])
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
[Disturbutionforce(1,[0,1],s->[0,-1e3+0*s,0],isgb)],
[Support(1,[true,true,false]),Support(2,[false,true,false])]
)
build_A(st2)
solvestruct(st2)
plot_orign(st2)
plot_moment(st2)
plot_shear(st2)