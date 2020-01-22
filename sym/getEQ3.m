function [ EQQQ ] = getEQ3(  )

EQQQ = sym(zeros(9,9,9));

parfor ij = 1:9
    [j,i] = ind2sub([3,3],ij);
    temp = sym(zeros(1,9,9));
    
    for k = 1:3
        for l = 1:3
            for m = 1:3
                for n = 1:3
                    temp(1,3*(k-1)+l,3*(m-1)+n) = EQ([i,j;k,l;m,n]);
                end
            end
        end
    end
    
    EQQQ(ij,:,:) = temp;
end

save('moments\EQQQSym','EQQQ');

end

