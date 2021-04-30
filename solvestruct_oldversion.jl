function heaviside(x::Float64)
    x*(x>0.0)
end

function nodeonele(i::Int64,els::Vector{Element})
    n = length(els)
    barindex = findall.([(els[ii].enode .== i) for ii = 1:n])
end

st = Struct(
    [Node([true,true,false],[0.0,0.0]),
    Node([false,true,false],[1.0,0.0]),
    Node([false,false,false],[0.0,1.0]),
    Node([false,false,false],[1.0,1.0])],
    [Element([1,2],1,1,1,[[true,true,true],[true,true,true]]),
    Element([1,3],1,1,1,[[true,true,true],[true,true,true]]),
    Element([2,3],1,1,1,[[true,true,true],[true,true,true]]),
    Element([3,4],1,1,1,[[true,true,true],[true,true,true]])],
    [Concentratedforce(4,0.5,[0,-1.0,0])],
)



node = st.globalnode;el = st.element;load = st.load

nvalue = 6*length(node)+12*length(el)


A = zeros(nvalue)
b = zeros(nvalue)

for i = 1:length(node)
    index = findall(vcat(node[i].state,.!(node[i].state)) )
    for j = 1:3
    temp = zeros(nvalue)
    temp[(i-1)*6 + index[j]] = 1
    A = hcat(A,temp)
    end
end

###
#A[1,1] = 1
#A[2,2] = 1
#A[3,6] = 1
#A[4,7] = 1
#A[5,11] = 1
#A[6,12] = 1

###connect
for i = 1:length(node)
    index = nodeonele(i,st.element)
    for j = 1:3
        index2 = findall((.!)(isempty.(index)))
        temp = zeros(nvalue,1+length(index2))
        
        if !(all(isempty.(index)))
             for (l,k) in enumerate(index2)
                if !isempty(index[k])
                    
                    temp[(k-1)*12+(index[k][1] - 1)*6+6*length(node)+j+3,1] = -1*Float64(el[k].constraint[(index[k][1])][j])
                    temp[(k-1)*12+(index[k][1] - 1)*6+6*length(node)+j,l+1] = -1*Float64(el[k].constraint[(index[k][1])][j])
                end
            end 
             
        end
        if !(iszero(temp))
            temp[(i-1)*6+j+3,1] = 1
            for l = 1:length(index2)
                temp[(i-1)*6+j,l+1] = 1
            end
            A = hcat(A,temp)
        else
            temp[(i-1)*6+j,1] = 1
            A = hcat(A,temp[:,1])
        end
    end
end
for i = 1:length(el)
    index = findall((.!)(hcat(el[i].constraint[1],el[i].constraint[2])))
    if !isempty(index)
        for j = 1:length(index)
        temp = zeros(nvalue)
        temp[(i-1)*12+(index[j][2]-1)*6+index[j][1]+3+6*length(node)] = 1
        A = hcat(A,temp)
        end
    end
end




#for i = 1:12
#A[6+i,i] = 1; A[6+i,i+12] = -1
#end


# equation of equilibrium
for i = 1:length(el)
    l = 1.0;
    for j = 1:3
        temp = zeros(nvalue)
        temp[(i-1)*12+6*length(node)+3+j] = 1
        temp[(i-1)*12+6*length(node)+9+j] = 1
        for k = 1:length(load)
            if load[k].elementcode == i
                b[size(A,2)] = -load[k].magnitude[j]
                if j == 3
                    temp[(i-1)*12+6*length(node)+9+2] = l
                    b[size(A,2)] = -load[k].magnitude[j] - load[k].magnitude[2]*load[k].loaction
                end
            end
        end
        A = hcat(A,temp)
    end
end


#A[19,16] = 1; A[19,18+4] = 1
#A[20,18] = 1; A[20,24] = 1;A[20,23]=1;b[20]=0.5;
#A[21,17] = 1; A[21,23] = 1;b[21] = 1

for i = 1:length(el)
    EA = el[i].E*el[i].A ;EI = el[i].E*el[i].I;
    node1 = el[i].enode[1];node2 = el[i].enode[2]
    gcn1 = st.globalnode[node1].coordinate;gcn2 = st.globalnode[node2].coordinate
    l = sqrt(transpose(gcn2-gcn1)*(gcn2-gcn1))
    ex_local = (gcn2-gcn1)./l
    ey_local = [-ex_local[2];ex_local[1]]
    

    for k = 1:length(load)
        F = load[k].magnitude[1:2]
        F_local = [transpose(ex_local);transpose(ey_local)]*F;
        loac = load[k].loaction
        if load[k].elementcode == i
            b[size(A,2)] = F[1]*(l-loac)
            b[size(A,2)+1] = -1/6*F[2]*(l-loac)^3 + 1/2*load[k].magnitude[3]*(l-loac)^2
            b[size(A,2)+2] = 1/2*F[2]*(l-loac)^2 - load[k].magnitude[3]*(l-loac)
        end
    end

    #u
    temp = zeros(nvalue)
    temp[(i-1)*12+6*length(node) .+ collect(1:2)] = ex_local
    temp[6+(i-1)*12+6*length(node) .+ collect(1:2)] = -ex_local
    temp[(i-1)*12+6*length(node) .+ collect(4:5)] = -l.*ex_local./EI
    A = hcat(A,temp) 
                
    #w
    temp = zeros(nvalue)
    temp[(i-1)*12+6*length(node) .+ collect(1:2)] = ey_local
    temp[(i-1)*12+6*length(node)+6 .+ collect(1:2)] = -ey_local
    temp[(i-1)*12+6*length(node)+3] = l
    temp[(i-1)*12+6*length(node)+6] = -1/2*l^2/EI
    temp[(i-1)*12+6*length(node) .+ collect(4:5)] = 1/6*l^2/EI*ey_local
    A = hcat(A,temp)
                

     #θ
    temp = zeros(nvalue)
    temp[(i-1)*12+6*length(node)+3] = 1
    temp[(i-1)*12+6*length(node)+6+3] = -1
    temp[(i-1)*12+6*length(node)+6] = l/EI   
    temp[(i-1)*12+6*length(node) .+ collect(4:5)] = -1/2*l^2/EI*ey_local
    A = hcat(A,temp)

end

row = size(A,1)
col = size(A,2)
B = zeros(row,col-1)
for i = 1:row
    for j = 2:col
        B[j-1,i] = A[i,j]
    end
end
            

#A[22,13] = 1;A[22,19] = -1;A[22,16] = -1
#A[23,14] = 1;A[23,20] = -1;A[23,17] = -1/6;A[23,15] = -1;b[23] = -1/6*0.5^3
#A[24,15] = 1;A[24,21] = -1;A[24,17] = 1/2;b[24] = 1/2*0.5^2
