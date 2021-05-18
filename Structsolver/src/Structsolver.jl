module Structsolver
export
    Element,Struct,Nodeforce,Concentratedforce,Disturbutionforce,Uniformforce,Support,Load
    build_A,build_frame,build_truss,
    forcecurves,internalforcedispsingle,leftsectionforcecurves,
    isgb,islc,
    plot_axis,plot_disp,plot_moment,plot_orign,plot_shear,
    solvestruct,
    totruss!

include("Structdefine.jl")
include("Structsolve.jl")
end # module
