A = zeros(24,24)
b = zeros(24)
A[1,1] = 1
A[2,2] = 1
A[3,6] = 1
A[4,7] = 1
A[5,8] = 1
A[6,12] = 1

for i = 1:12
A[6+i,i] = 1; A[6+i,i+12] = -1
end

A[19,16] = 1; A[19,18+4] = 1
A[20,18] = 1; A[20,24] = 1;A[20,23]=1;b[20]=0.5;
A[21,17] = 1; A[21,23] = 1;b[21] = 1

A[22,13] = 1;A[22,19] = -1;A[22,16] = -1
A[23,14] = 1;A[23,20] = -1;A[23,17] = -1/6;A[23,15] = -1;b[23] = -1/6*0.5^3
A[24,15] = 1;A[24,21] = -1;A[24,17] = 1/2;b[24] = 1/2*0.5^2


#buck
for k = 1:length(index)
    #@show el[k].constraint[index[k][1]][j]
    #@show typeof(el[k].constraint[index[k][1]][j])
    if !isempty(index[k])
        if (el[k].constraint[index[k][1]][j]) #displacement is same
        temp = zeros(nvalue)
        temp[sum(svalue)+(k-1)*12+j] = 1
        A = hcat(A,temp)
        end
    end
end