function [ output ] = dsdT0( s, inds, indT )

U = eye(3);
V = eye(3);

if size(indT,1)==1
    output = U(indT(1,1),inds)*V(indT(1,2),inds);
else
    N = size(indT(2:end,:),1);
    output = U(indT(1,1),inds)*dVdT0(s,[indT(1,2),inds],indT(2:end,:));
    
    for n = 1:N
        comb = combnk(1:N,n);
        Nc = size(comb,1);
        
        for nc = 1:Nc
            dUInd = comb(nc,:)+1;
            dVInd = setdiff(1:N,comb(nc,:))+1;
            
            if n==N
                output = output + dUdT0(s,[indT(1,1),inds],indT(2:end,:))*V(indT(1,2),inds);
            else
                output = output + dUdT0(s,[indT(1,1),inds],indT(dUInd,:))...
                    *dVdT0(s,[indT(1,2),inds],indT(dVInd,:));
            end
        end
    end
end

end


function [ output ] = dUdT0( s, indU, indT )

U = eye(3);
V = eye(3);

ext = setdiff([1,2,3],indU(2));

i = indT(1,1);
j = indT(1,2);
k = ext(1);
l = indU(2);
OmegaU1 = (s(l)*U(i,k)*V(j,l)+s(k)*U(i,l)*V(j,k))/(s(l)^2-s(k)^2);

k = ext(2);
OmegaU2 = (s(l)*U(i,k)*V(j,l)+s(k)*U(i,l)*V(j,k))/(s(l)^2-s(k)^2);

if size(indT,1)==1
    output = U(indU(1),ext(1))*OmegaU1 + U(indU(1),ext(2))*OmegaU2;
else
    N = size(indT(2:end,:),1);
    
    % DUOmegaU1
    DUOmegaU1 = U(indU(1),ext(1))*dOmegaUdT0(s,indT(1,:),[ext(1),indU(2)],indT(2:end,:));
    for n = 1:N
        comb = combnk(1:N,n);
        Nc = size(comb,1);
        
        for nc = 1:Nc
            dUInd = comb(nc,:)+1;
            dOmegaUInd = setdiff(1:N,comb(nc,:))+1;
            
            if n==N
                DUOmegaU1 = DUOmegaU1 + dUdT0(s,[indU(1),ext(1)],indT(2:end,:))*OmegaU1;
            else
                DUOmegaU1 = DUOmegaU1 + dUdT0(s,[indU(1),ext(1)],indT(dUInd,:))...
                    *dOmegaUdT0(s,indT(1,:),[ext(1),indU(2)],indT(dOmegaUInd,:));
            end
        end
    end
    
    % DUOmegaU2
    DUOmegaU2 = U(indU(1),ext(2))*dOmegaUdT0(s,indT(1,:),[ext(2),indU(2)],indT(2:end,:));
    for n = 1:N
        comb = combnk(1:N,n);
        Nc = size(comb,1);
        
        for nc = 1:Nc
            dUInd = comb(nc,:)+1;
            dOmegaUInd = setdiff(1:N,comb(nc,:))+1;
            
            if n==N
                DUOmegaU2 = DUOmegaU2 + dUdT0(s,[indU(1),ext(2)],indT(2:end,:))*OmegaU2;
            else
                DUOmegaU2 = DUOmegaU2 + dUdT0(s,[indU(1),ext(2)],indT(dUInd,:))...
                    *dOmegaUdT0(s,indT(1,:),[ext(2),indU(2)],indT(dOmegaUInd,:));
            end
        end
    end
    
    output = DUOmegaU1+DUOmegaU2;
end

end


function [ output ] = dVdT0( s, indV, indT )

U = eye(3);
V = eye(3);

ext = setdiff([1,2,3],indV(2));

i = indT(1,1);
j = indT(1,2);
k = ext(1);
l = indV(2);
OmegaV1 = (s(k)*U(i,k)*V(j,l)+s(l)*U(i,l)*V(j,k))/(s(k)^2-s(l)^2);

k = ext(2);
OmegaV2 = (s(k)*U(i,k)*V(j,l)+s(l)*U(i,l)*V(j,k))/(s(k)^2-s(l)^2);

if size(indT,1)==1
    output = -V(indV(1),ext(1))*OmegaV1 - V(indV(1),ext(2))*OmegaV2;
else
    N = size(indT(2:end,:),1);
    
    % DVOmegaV1
    DVOmegaV1 = V(indV(1),ext(1))*dOmegaVdT0(s,indT(1,:),[ext(1),indV(2)],indT(2:end,:));
    for n = 1:N
        comb = combnk(1:N,n);
        Nc = size(comb,1);
        
        for nc = 1:Nc
            dVInd = comb(nc,:)+1;
            dOmegaVInd = setdiff(1:N,comb(nc,:))+1;
            
            if n==N
                DVOmegaV1 = DVOmegaV1 + dVdT0(s,[indV(1),ext(1)],indT(2:end,:))*OmegaV1;
            else
                DVOmegaV1 = DVOmegaV1 + dVdT0(s,[indV(1),ext(1)],indT(dVInd,:))...
                    *dOmegaVdT0(s,indT(1,:),[ext(1),indV(2)],indT(dOmegaVInd,:));
            end
        end
    end
    
    % DUOmegaU2
    DVOmegaV2 = V(indV(1),ext(2))*dOmegaVdT0(s,indT(1,:),[ext(2),indV(2)],indT(2:end,:));
    for n = 1:N
        comb = combnk(1:N,n);
        Nc = size(comb,1);
        
        for nc = 1:Nc
            dVInd = comb(nc,:)+1;
            dOmegaVInd = setdiff(1:N,comb(nc,:))+1;
            
            if n==N
                DVOmegaV2 = DVOmegaV2 + dVdT0(s,[indV(1),ext(2)],indT(2:end,:))*OmegaV2;
            else
                DVOmegaV2 = DVOmegaV2 + dVdT0(s,[indV(1),ext(2)],indT(dVInd,:))...
                    *dOmegaVdT0(s,indT(1,:),[ext(2),indV(2)],indT(dOmegaVInd,:));
            end
        end
    end
    
    output = -DVOmegaV1-DVOmegaV2;
end

end


function [ output ] = dOmegaUdT0( s, indOsup, indOsub, indT )

U = eye(3);
V = eye(3);

i = indOsup(1);
j = indOsup(2);
k = indOsub(1);
l = indOsub(2);
OmegaU = (s(l)*U(i,k)*V(j,l)+s(k)*U(i,l)*V(j,k))/(s(l)^2-s(k)^2);
OmegaV = (s(k)*U(i,k)*V(j,l)+s(l)*U(i,l)*V(j,k))/(s(k)^2-s(l)^2);

N = size(indT,1);

% DslOmegaUk1
DslOmegaUkl = 0;
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    for nc = 1:Nc
        dsInd = comb(nc,:);
        dOmegaInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DslOmegaUkl = DslOmegaUkl + dsdT0(s,indOsub(2),indT)*OmegaU;
        else
            DslOmegaUkl = DslOmegaUkl + dsdT0(s,indOsub(2),indT(dsInd,:))...
                *dOmegaUdT0(s,indOsup,indOsub,indT(dOmegaInd,:));
        end
    end
end

% DskOmegaVkl
DskOmegaVkl = 0;
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    for nc = 1:Nc
        dsInd = comb(nc,:);
        dOmegaInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DskOmegaVkl = DskOmegaVkl + dsdT0(s,indOsub(1),indT)*OmegaV;
        else
            DskOmegaVkl = DskOmegaVkl + dsdT0(s,indOsub(1),indT(dsInd,:))...
                *dOmegaVdT0(s,indOsup,indOsub,indT(dOmegaInd,:));
        end
    end
end

%DskOmegaUkl
DskOmegaUkl = 0;
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    for nc = 1:Nc
        dsInd = comb(nc,:);
        dOmegaInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DskOmegaUkl = DskOmegaUkl + dsdT0(s,indOsub(1),indT)*OmegaU;
        else
            DskOmegaUkl = DskOmegaUkl + dsdT0(s,indOsub(1),indT(dsInd,:))...
                *dOmegaUdT0(s,indOsup,indOsub,indT(dOmegaInd,:));
        end
    end
end

% DslOmegaVkl
DslOmegaVkl = 0;
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    for nc = 1:Nc
        dsInd = comb(nc,:);
        dOmegaInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DslOmegaVkl = DslOmegaVkl + dsdT0(s,indOsub(2),indT)*OmegaV;
        else
            DslOmegaVkl = DslOmegaVkl + dsdT0(s,indOsub(2),indT(dsInd,:))...
                *dOmegaVdT0(s,indOsup,indOsub,indT(dOmegaInd,:));
        end
    end
end

% DUikVjl
DUikVjl = U(indOsup(1),indOsub(1))*dVdT0(s,[indOsup(2),indOsub(2)],indT);
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    
    for nc = 1:Nc
        dUInd = comb(nc,:);
        dVInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DUikVjl = DUikVjl + dUdT0(s,[indOsup(1),indOsub(1)],indT)*...
                V(indOsup(2),indOsub(2));
        else
            DUikVjl = DUikVjl + dUdT0(s,[indOsup(1),indOsub(1)],indT(dUInd,:))*...
                dVdT0(s,[indOsup(2),indOsub(2)],indT(dVInd,:));
        end
    end
end

% DUilVjk
DUilVjk = U(indOsup(1),indOsub(2))*dVdT0(s,[indOsup(2),indOsub(1)],indT);
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    
    for nc = 1:Nc
        dUInd = comb(nc,:);
        dVInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DUilVjk = DUilVjk + dUdT0(s,[indOsup(1),indOsub(2)],indT)*...
                V(indOsup(2),indOsub(1));
        else
            DUilVjk = DUilVjk + dUdT0(s,[indOsup(1),indOsub(2)],indT(dUInd,:))*...
                dVdT0(s,[indOsup(2),indOsub(1)],indT(dVInd,:));
        end
    end
end

% solve the equation
A = [s(indOsub(2)),s(indOsub(1));s(indOsub(1)),s(indOsub(2))];
y = [DUikVjl-DslOmegaUkl-DskOmegaVkl;
    -DUilVjk-DskOmegaUkl-DslOmegaVkl];

DOmega = linsolve(A,y);
output = DOmega(1);

end


function [ output ] = dOmegaVdT0( s, indOsup, indOsub, indT )

U = eye(3);
V = eye(3);

i = indOsup(1);
j = indOsup(2);
k = indOsub(1);
l = indOsub(2);
OmegaU = (s(l)*U(i,k)*V(j,l)+s(k)*U(i,l)*V(j,k))/(s(l)^2-s(k)^2);
OmegaV = (s(k)*U(i,k)*V(j,l)+s(l)*U(i,l)*V(j,k))/(s(k)^2-s(l)^2);

N = size(indT,1);

% DslOmegaUk1
DslOmegaUkl = 0;
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    for nc = 1:Nc
        dsInd = comb(nc,:);
        dOmegaInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DslOmegaUkl = DslOmegaUkl + dsdT0(s,indOsub(2),indT)*OmegaU;
        else
            DslOmegaUkl = DslOmegaUkl + dsdT0(s,indOsub(2),indT(dsInd,:))...
                *dOmegaUdT0(s,indOsup,indOsub,indT(dOmegaInd,:));
        end
    end
end

% DskOmegaVkl
DskOmegaVkl = 0;
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    for nc = 1:Nc
        dsInd = comb(nc,:);
        dOmegaInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DskOmegaVkl = DskOmegaVkl + dsdT0(s,indOsub(1),indT)*OmegaV;
        else
            DskOmegaVkl = DskOmegaVkl + dsdT0(s,indOsub(1),indT(dsInd,:))...
                *dOmegaVdT0(s,indOsup,indOsub,indT(dOmegaInd,:));
        end
    end
end

%DskOmegaUkl
DskOmegaUkl = 0;
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    for nc = 1:Nc
        dsInd = comb(nc,:);
        dOmegaInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DskOmegaUkl = DskOmegaUkl + dsdT0(s,indOsub(1),indT)*OmegaU;
        else
            DskOmegaUkl = DskOmegaUkl + dsdT0(s,indOsub(1),indT(dsInd,:))...
                *dOmegaUdT0(s,indOsup,indOsub,indT(dOmegaInd,:));
        end
    end
end

% DslOmegaVkl
DslOmegaVkl = 0;
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    for nc = 1:Nc
        dsInd = comb(nc,:);
        dOmegaInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DslOmegaVkl = DslOmegaVkl + dsdT0(s,indOsub(2),indT)*OmegaV;
        else
            DslOmegaVkl = DslOmegaVkl + dsdT0(s,indOsub(2),indT(dsInd,:))...
                *dOmegaVdT0(s,indOsup,indOsub,indT(dOmegaInd,:));
        end
    end
end

% DUikVjl
DUikVjl = U(indOsup(1),indOsub(1))*dVdT0(s,[indOsup(2),indOsub(2)],indT);
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    
    for nc = 1:Nc
        dUInd = comb(nc,:);
        dVInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DUikVjl = DUikVjl + dUdT0(s,[indOsup(1),indOsub(1)],indT)*...
                V(indOsup(2),indOsub(2));
        else
            DUikVjl = DUikVjl + dUdT0(s,[indOsup(1),indOsub(1)],indT(dUInd,:))*...
                dVdT0(s,[indOsup(2),indOsub(2)],indT(dVInd,:));
        end
    end
end

% DUilVjk
DUilVjk = U(indOsup(1),indOsub(2))*dVdT0(s,[indOsup(2),indOsub(1)],indT);
for n = 1:N
    comb = combnk(1:N,n);
    Nc = size(comb,1);
    
    for nc = 1:Nc
        dUInd = comb(nc,:);
        dVInd = setdiff(1:N,comb(nc,:));
        
        if n==N
            DUilVjk = DUilVjk + dUdT0(s,[indOsup(1),indOsub(2)],indT)*...
                V(indOsup(2),indOsub(1));
        else
            DUilVjk = DUilVjk + dUdT0(s,[indOsup(1),indOsub(2)],indT(dUInd,:))*...
                dVdT0(s,[indOsup(2),indOsub(1)],indT(dVInd,:));
        end
    end
end

% solve the equation
A = [s(indOsub(2)),s(indOsub(1));s(indOsub(1)),s(indOsub(2))];
y = [DUikVjl-DslOmegaUkl-DskOmegaVkl;
    -DUilVjk-DskOmegaUkl-DslOmegaVkl];

DOmega = linsolve(A,y);
output = DOmega(2);

end

