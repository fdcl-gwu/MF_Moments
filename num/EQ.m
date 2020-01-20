function [ output ] = EQ( s, ind )

n = size(ind,1);

output = 0;
for k = 1:n
    comb = setC(n,k);
    Nc = length(comb);
    
    for nc = 1:Nc
        % detect if there is an index set l such that \partial_l s vanishes
        isZero = cellfun(@(x)sum(rem(histcounts(ind(x,:)),2)),comb{nc});
        if ~isempty(find(isZero,1))
            continue;
        end
        
        DcDsInd = {};
        DsDT = [];
        for nl = 1:k
            [DcDsInd,DsDT] = sumAlpha(s,DcDsInd,DsDT,ind(comb{nc}{nl},:));
        end
        
        DcDs = zeros(3^k,1);
        for nk = 1:3^k
            if DsDT(nk)==0
                % if DsDT=0, no need to calculate DcDs
                % This is crucial to calculate dcds recursively! Otherwise
                % the recursion will be infinite!
                DcDs(nk) = 0;
            else
                DcDs(nk) = dcds(s,DcDsInd{nk});
            end
        end
        
        output = output+sum(DcDs.*DsDT);
    end
end

end


function [ comb ] = setC( n, k )

if n==1 && k==1
    comb = {{1}};
    return;
elseif n==2 && k==1
    comb = {{[1,2]}};
    return;
end

if k>1
    comb = setC(n-1,k-1);
    Nc = length(comb);
    for nc = 1:Nc
        comb{nc} = [comb{nc},n];
    end
end

if n>k
    comb2 = setC(n-1,k);
    Nc = length(comb2);
    Nk = k;
    for nc = 1:Nc
        for nk = 1:Nk
            tempComb = comb2{nc};
            tempComb{nk} = [tempComb{nk},n];
            if ~exist('comb','var')
                comb = {tempComb};
            else
                comb = [comb,{tempComb}]; %#ok<AGROW>
            end
        end
    end
end

end


function [ DcDsInd, DsDT ] = sumAlpha( s, DcDsInd, DsDT, l )

if isempty(DcDsInd)
    for nAlpha = 1:3
        DcDsInd{nAlpha,1} = nAlpha;
        DsDT(nAlpha,1) = dsdT0(s,nAlpha,l);
    end
else
    N = length(DcDsInd);
    tempDcDsInd = DcDsInd;
    tempDsDT = DsDT;
    DcDsInd = cell(N*3,1);
    DsDT = zeros(N*3,1);
    
    for n = 1:N
        for nAlpha = 1:3
            DcDsInd{3*(n-1)+nAlpha} = [tempDcDsInd{n};nAlpha];
            DsDT(3*(n-1)+nAlpha) = tempDsDT(n)*dsdT0(s,nAlpha,l);
        end
    end
end

end

