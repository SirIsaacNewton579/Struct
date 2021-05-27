module Structsolver
export
    Element,Struct,Nodeforce,Concentratedforce,Disturbutionforce,Uniformforce,Support,Load,
    build_A,build_frame,build_truss,totruss!,build_fix,build_hinge
    forcecurves,internalforcedispsingle,leftsectionforcecurves,
    isgb,islc,
    plot_axis,plot_disp,plot_moment,plot_orign,plot_shear,
    solvestruct,
    
include("Structdefine.jl")
include("Structsolve.jl")
end # module
