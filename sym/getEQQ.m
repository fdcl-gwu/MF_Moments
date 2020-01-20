function [ EQQ ] = getEQQ(  )

EQQ = sym(zeros(9,9));

parfor ij = 1:9
    [j,i] = ind2sub([3,3],ij);
    temp = sym(zeros(1,9));
    for k = 1:3
        for l = 1:3
            temp(3*(k-1)+l) = EQ([i,j;k,l]);
        end
    end
    EQQ(ij,:) = temp;
end

save('moments\EQQSym','EQQ');

end

