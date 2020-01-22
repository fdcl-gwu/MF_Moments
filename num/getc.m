function [ c ] = getc( s )

c = integral(@(u) f_kunze_s(u,s),-1,1);

end


function Y=f_kunze_s(u,s)
% integrand for the normalizing constant
    J=besseli(0,1/2*(s(1)-s(2))*(1-u)).*besseli(0,1/2*(s(1)+s(2))*(1+u));
    Y=1/2*exp(s(3)*u).*J;
end

