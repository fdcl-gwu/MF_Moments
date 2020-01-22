function [ output ] = dcds( s, ind )

if length(ind)==1
    output = integral(@(u) f_kunze_s_deriv_i(u,s,ind),-1,1);
elseif length(ind)==2
    if ind(1)==ind(2)
        c = getc(s);
        
        ext = setdiff([1,2,3],ind);
        indj = [ind(1),ext(1);ind(1),ext(1)];
        indk = [ind(1),ext(2);ind(1),ext(2)];
        output = c-EQ(s,indj)-EQ(s,indk);
    else
        ext = setdiff([1,2,3],ind);
        output = EQ(s,[ext,ext])+EQ(s,[ind';flip(ind')]);
    end
else
    ind1 = repmat(ind(3:end),1,2);
    if ind(1)==ind(2)
        ext = setdiff([1,2,3],ind(1:2));
        indj = [ind(1),ext(1);ind(1),ext(1);ind1];
        indk = [ind(1),ext(2);ind(1),ext(2);ind1];
        output = EQ(s,ind1) - EQ(s,indj) - EQ(s,indk);
    else
        ext = setdiff([1,2,3],ind(1:2));
        output = EQ(s,[ext,ext;ind1])+EQ(s,[ind(1:2)';flip(ind(1:2)');ind1]);
    end
end

end


function Y=f_kunze_s_deriv_i(u,s,i)
% integrand for the derivative of the normalizing constant
index=circshift([1 2 3],[0 4-i]);
j=index(2);
k=index(3);

J00=besseli(0,1/2*(s(j)-s(k))*(1-u)).*besseli(0,1/2*(s(j)+s(k))*(1+u));
Y=1/2*J00.*u.*exp(s(i)*u);
end

