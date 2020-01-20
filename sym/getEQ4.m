function [ EQ4 ] = getEQ4(  )

EQ4 = sym(zeros(9,9,9,9));

parfor ij1 = 1:9
    [j1,i1] = ind2sub([3,3],ij1);
    temp = sym(zeros(1,9,9,9));
    
    for i2 = 1:3
        for j2 = 1:3
            for i3 = 1:3
                for j3 = 1:3
                    for i4 = 1:3
                        for j4 = 1:3
                            temp(1,3*(i2-1)+j2,3*(i3-1)+j3,3*(i4-1)+j4) = EQ([i1,j1;i2,j2;i3,j3;i4,j4]);
                        end
                    end
                end
            end
        end
    end
    
    EQ4(ij1,:,:,:) = temp;
end

save('\moments\EQ4Sym','EQ4');

end

