using TaylorSeries ; Ts = TaylorSeries
using SymPy; Sp = SymPy
function nodeonele(i::Int64,els::Vector{Element})
    n = length(els)
    barindex = findall.([(els[ii].enode .== i) for ii = 1:n])
end

function forcecurves(F :: Function,x1,x2,l; Fₗ=[0.0;0.0;0.0])
Fx = Fₗ[1];Fy = Fₗ[2];M = Fₗ[3];
    try 
        t = Taylor1(Float64, 10)
        p0 = F(t) 
        p = Vector{Vector{Taylor1{Float64}}}(undef,4)
        pnapd = Vector{Vector{Taylor1{Float64}}}(undef,4)
        ∫F = Vector{Vector{Taylor1{Float64}}}(undef,4)
        
        p[1] = Ts.integrate.(p0)+[Fx;Fy;0]
        ∫F[1] = p[1]-p[1](x1)
        pnapd[1] = [0,0,0]*t + ∫F[1](x2)
        for jj = 2:4
        p[jj] = Ts.integrate.(∫F[jj-1])
        (jj==2 ? p[jj][3] += -M : nothing)
        ∫F[jj] = p[jj] - p[jj](x1)
        pnapd[jj] = Ts.integrate.(pnapd[jj-1]) -Ts.integrate.(pnapd[jj-1])(x2) + ∫F[jj](x2)
        end
        [Fx;Fy;-M+Fy*t],∫F,pnapd 
    catch
        @warn "can't obtain TaylorSeries at 0.0,try analytic solution."
        try
            x = Sp.symbols("x")
            i1p1 = Sp.integrate.(F,x1,x)+[Fx;Fy;0]; i1p2 = Sp.integrate.(F,x2,x); p1apds = i1p1 - i1p2;
            i2p1 = Sp.integrate.(i1p1,x1,x)-M; i2p2 = Sp.integrate.(i1p2,x2,x); p2apds = i2p1 - i2p2;
            i3p1 = Sp.integrate.(i2p1,x1,x); i3p2 = Sp.integrate.(i2p2,x2,x); p3apds = i3p1 - i3p2;
            i4p1 = Sp.integrate.(i3p1,x1,x); i4p2 = Sp.integrate.(i3p2,x2,x); p4apds = i4p1 - i4p2;
            inp1 = [i1p1,i2p1,i3p1,i4p1]
            pnapds = [p1apds,p2apds,p3apds,p4apds]
            ∫F = Vector{Function}(undef,4)
            pnapd = Vector{Function}(undef,4)
            for m = 1:4 
                hfF = [lambdify(inp1[m][i]) for i = 0:2]
                ∫F[m] = x->[hfF[i](x) for i = 1:3]
                hfpa = [lambdify(pnapds[m][i]) for i = 0:2]
                pnapd = x->[hfpa[i](x) for i = 1:3]
            end
            [Fx;Fy;-M+Fy*t],∫F,pnapd
        catch
            t = Taylor1(Float64, 10)
            p1 = taylor_expand(F,x1) 
            p2 = taylor_expand(F,x2) 
            i1p1 = Ts.integrate.(p1); i1p2 = Ts.integrate.(p2); p1apd = i1p1(t-x1) - i1p2(t-x2);
            i2p1 = Ts.integrate.(i1p1); i2p2 = Ts.integrate.(i1p2); p2apd = i2p1(t-x1) - i2p2(t-x2);
            i3p1 = Ts.integrate.(i2p1); i3p2 = Ts.integrate.(i2p2); p3apd = i3p1(t-x1) - i3p2(t-x2);
            i4p1 = Ts.integrate.(i3p1); i4p2 = Ts.integrate.(i3p2); p4apd = i4p1(t-x1) - i4p2(t-x2);
            pnapd =  [p1apd,p2apd,p3apd,p4apd]
            nothing,nothing,pnapd 
        end
    end
end

function build_A_consider_support!(A,sp :: Vector{Support}, el :: Vector{Element})
    spb = [sup.ncst for sup in sp]
    spbindex = findall.(spb)
    svalue = length.(spbindex)
    nvalue = 12*length(el) + sum(svalue)
    for i = 1:length(sp)
        index = nodeonele(sp[i].snode,el)
        if !(all(isempty.(index)))
            for j in spbindex[i]
                k = findall((.!)(isempty.(index)))
                n = findfirst([el[m].constraint[index[m][1]][j] for m in k])
                temp = zeros(nvalue+1)
                temp[sum(svalue)+(k[n]-1)*12+j+(index[k[n]][1]-1)*6] = 1  #displacement
                push!(A,temp)
            end
        end
    end

    for i = 1:length(sp)
        index = nodeonele(sp[i].snode,el)
        if !(all(isempty.(index)))
            for j = 1:3
                temp = zeros(nvalue+1)
                ord = 0
                if spb[i][j]
                    if i == 1
                        
                        #@show [spbindex[i] .== j]
                        ord = findfirst(spbindex[i] .== j)
                    else
                        ord = sum(svalue[1:i-1]) + findfirst(spbindex[i] .== j)
                    end
                    #@show ord
                    temp[ord] = 1
                end
                for k = 1:length(index)
                    if !isempty(index[k])
                        if (el[k].constraint[index[k][1]][j])
                            temp[sum(svalue)+(k-1)*12+3+j+(index[k][1]-1)*6] = -1   #force
                        end
                    end
                end
                #@show temp
                (!iszero(temp) ? push!(A,temp) : nothing)
            end
        end
    end
end


function build_A_consider_connect!(A,node :: Vector{Vector{Float64}},el :: Vector{Element},sp :: Vector{Support},nvalue,svalue)
    #connect
    for i = 1:length(node)
        index = nodeonele(i,el)
        k = findall((.!)(isempty.(index)))
        for j = 1:3
            m = findall([el[m].constraint[index[m][1]][j] for m in k])
            if length(m)>1
                for n in m[2:end]
                    temp = zeros(nvalue+1)
                    temp[sum(svalue)+(k[m[1]]-1)*12+j+(index[k[m[1]]][1]-1)*6] = 1
                    temp[sum(svalue)+(k[n]-1)*12+j+(index[k[n]][1]-1)*6] = -1
                    push!(A,temp)
                end
            end
        end
        
        if !any([spi.snode == i for spi in sp])
            for j = 1:3
                m = findall([el[m].constraint[index[m][1]][j] for m in k])
                temp = zeros(nvalue+1)
                for n in m 
                    temp[sum(svalue)+(k[n]-1)*12+j+3+(index[k[n]][1]-1)*6] = 1  #force
                end
                push!(A,temp)
            end
        end
    end

    for i = 1:length(el)
        index = findall((.!)(hcat(el[i].constraint[1],el[i].constraint[2])))
        if !isempty(index)
            for j = 1:length(index)
                temp = zeros(nvalue+1)
                temp[sum(svalue) + (i-1)*12+(index[j][2]-1)*6+index[j][1]+3] = 1  #force
                push!(A,temp)
            end
        end
    end
end



function single_element_equilibrium!(A,
    node :: Vector{Vector{Float64}} ,
    eli :: Element ,
    i :: Int64,
    load :: Vector{Load},
    nvalue,svalue)

    EA = eli.E*eli.A ;EI = eli.E*eli.I;
    node1 = eli.enode[1];node2 = eli.enode[2]
    gcn1 = node[node1];gcn2 = node[node2]
    l = norm(gcn2-gcn1)
    ex_local = (gcn2-gcn1)/l
    ey_local = [-ex_local[2];ex_local[1]]
    T = [transpose(ex_local) 0;transpose(ey_local) 0; 0 0 1]
           
    function consider_load_equilibrium!(loadi :: Concentratedforce)
        F = loadi.fgol ? loadi.magnitude : T\loadi.magnitude
        F_local = T*F;
        loac = loadi.loaction
        b[1:3] +=  -F
        b[3] +=  F_local[2]*loac*l
        
    end
    function consider_load_equilibrium!(loadi :: Disturbutionforce)
        loac = loadi.loaction
        x1 =  loac[1]*l;x2 = loac[2]*l;
        #@show l
        F = loadi.fgol ? loadi.magnitude : x->T\loadi.magnitude(x)
        
        ~,~,pnapd = forcecurves(F,x1,x2,l)
        
        
        pnapd =  [(i(l)) for i in pnapd]
        num_∫ⁿF =  [T*i for i in pnapd]
        
        b[1:3] +=  -pnapd[1]
        b[3] +=  num_∫ⁿF[2][2]
        @show b
        
    end
    function consider_load_equilibrium!(loadi :: Uniformforce)
        loac = loadi.loaction
        x1 =  loac[1]*l;x2 = loac[2]*l;
        F = loadi.fgol ? loadi.magnitude : T\loadi.magnitude
        ∫F = F*(x2-x1)
        ∫²F = F/2*((l-x1)^2-(l-x2)^2)
        ∫³F = F/6*((l-x1)^3-(l-x2)^3)
        ∫⁴F = F/24*((l-x1)^4-(l-x2)^4)
        
        num_∫F = T*∫F
        num_∫²F = T*∫²F
        num_∫³F = T*∫³F
        num_∫⁴F = T*∫⁴F
        
        b[1:3] +=  -∫F
        b[3] +=  num_∫²F[2]
        
    end

    b = zeros(3)
    for k = 1:length(load)
        if load[k].elementcode == i
        consider_load_equilibrium!(load[k])
        end
    end
    
    for j = 1:3
        temp = zeros(nvalue+1)
        temp[(i-1)*12+sum(svalue)+3+j] = 1
        temp[(i-1)*12+sum(svalue)+9+j] = 1
        if j == 3
            temp[(i-1)*12+sum(svalue)+3 .+ collect(1:2)] = -l*ey_local
        end
        temp[end] = b[j]
        push!(A,temp)
    end
end

function single_element_displacement!(A,
    node :: Vector{Vector{Float64}} ,
    eli :: Element ,
    i :: Int64,
    load :: Vector{Load},
    nvalue,svalue)

    EA = eli.E*eli.A ;EI = eli.E*eli.I;
    node1 = eli.enode[1];node2 = eli.enode[2]
    gcn1 = node[node1];gcn2 = node[node2]
    l = norm(gcn2-gcn1)
    ex_local = (gcn2-gcn1)/l
    ey_local = [-ex_local[2];ex_local[1]]
    T = [transpose(ex_local) 0;transpose(ey_local) 0; 0 0 1]
    function consider_load_displacement!(loadi :: Concentratedforce)
        F = loadi.fgol ? loadi.magnitude : T\loadi.magnitude
        F_local = T*F;
        loac = loadi.loaction
        b[1] = F_local[1]*(l-loac*l)/EA
        b[2] = (-1/6*F_local[2]*(l-loac*l)^3 + 1/2*loadi.magnitude[3]*(l-loac*l)^2)/EI
        b[3] = (-1/2*F_local[2]*(l-loac*l)^2 - loadi.magnitude[3]*(l-loac*l))/EI
    end
    function consider_load_displacement!(loadi :: Disturbutionforce)
        loac = loadi.loaction
        x1 =  loac[1]*l;x2 = loac[2]*l;
        F = loadi.fgol ? loadi.magnitude : x->T\loadi.magnitude(x)

        ~,~,pnapd = forcecurves(F,x1,x2,l)

        pnapd =  [(i(l)) for i in pnapd]
        num_∫ⁿF =  [T*i for i in pnapd]
        b[1] +=  num_∫ⁿF[2][1]/EA  #u
        
        b[2] += (-num_∫ⁿF[4][2] + num_∫ⁿF[3][3])/EI #v
        b[3] +=  (-num_∫ⁿF[3][2] + num_∫ⁿF[2][3])/EI  #θ
        # @show b
    end
    function consider_load_displacement!(loadi :: Uniformforce)
        loac = loadi.loaction
        x1 =  loac[1]*l;x2 = loac[2]*l;
        F = loadi.fgol ? loadi.magnitude : T\loadi.magnitude
        ∫F = F*(x2-x1)
        ∫²F = F/2*((l-x1)^2-(l-x2)^2)
        ∫³F = F/6*((l-x1)^3-(l-x2)^3)
        ∫⁴F = F/24*((l-x1)^4-(l-x2)^4)

        num_∫F = T*∫F
        num_∫²F = T*∫²F
        num_∫³F = T*∫³F
        num_∫⁴F = T*∫⁴F
        
        b[1] +=  num_∫²F[1]/EA  #u
        
        b[2] += (-num_∫⁴F[2] + num_∫³F[3])/EI #v
        b[3] +=  (-num_∫³F[2] + num_∫²F[3])/EI  #θ
    end

    b = zeros(3,1)
    for k = 1:length(load)
        if load[k].elementcode == i
            consider_load_displacement!(load[k])
        end
    end

    #u
    temp = zeros(nvalue+1)
    temp[(i-1)*12+sum(svalue) .+ collect(1:2)] = ex_local
    temp[6+(i-1)*12+sum(svalue) .+ collect(1:2)] = -ex_local
    temp[(i-1)*12+sum(svalue) .+ collect(4:5)] = -l.*ex_local./EA
    temp[end]=b[1]
    push!(A,temp) 
                
    #w
    temp = zeros(nvalue+1)
    temp[(i-1)*12+sum(svalue) .+ collect(1:2)] = ey_local
    temp[(i-1)*12+sum(svalue)+6 .+ collect(1:2)] = -ey_local
    temp[(i-1)*12+sum(svalue)+3] = l
    temp[(i-1)*12+sum(svalue)+6] = -1/2*l^2/EI
    temp[(i-1)*12+sum(svalue) .+ collect(4:5)] = 1/6*l^3/EI*ey_local
    temp[end]=b[2]
    push!(A,temp)
                

    #θ
    temp = zeros(nvalue+1)
    temp[(i-1)*12+sum(svalue)+3] = 1
    temp[(i-1)*12+sum(svalue)+6+3] = -1
    temp[(i-1)*12+sum(svalue)+6] = -l/EI   
    temp[(i-1)*12+sum(svalue) .+ collect(4:5)] = +1/2*l^2/EI*ey_local
    temp[end]=b[3]
    push!(A,temp)

end

function build_A_consider_equilibrium_and_displacement!(A,
    node :: Vector{Vector{Float64}} ,
    el :: Vector{Element} ,
    load :: Vector{Load},
    nvalue,svalue)
    # equilibrium equations and displacement equations
    for i = 1:length(el)
        
        single_element_equilibrium!(A,node,el[i],i,load,nvalue,svalue)

        single_element_displacement!(A,node,el[i],i,load,nvalue,svalue)
    end
end