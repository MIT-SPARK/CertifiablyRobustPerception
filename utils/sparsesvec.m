function At   = sparsesvec(blk,Acell)
%% Implement svec of a cell of sparse matrices
%% This is much faster than calling svec(blk,Acell) using mexsvec

n           = blk{1,2};
m           = length(Acell);
ndelta      = n*(n+1)/2;
rsq2        = sqrt(2);

ssi         = [];
ssj         = [];
ssdata      = [];
% fprintf('svec a cell of sparse matrices, progress ...')
for i = 1:m
    % if rem(i,10000) == 1
    %     fprintf('%d/%d ',i,m);
    % end
    Ai              = Acell{i};
    [ii,jj,val]     = find(Ai);
    mask_triu       = (jj >= ii);
    ii              = ii(mask_triu);
    jj              = jj(mask_triu);
    val             = val(mask_triu);
    
    mask_diag       = (ii == jj);
    mask_nondiag    = ~mask_diag;
    
    si              = ((jj-1).*jj)./2 + ii;
    sj              = ones(length(si),1) * i;
    sdata           = val;
    sdata(mask_nondiag) = sdata(mask_nondiag) * rsq2;
    
    ssi             = [ssi;si(:)];
    ssj             = [ssj;sj(:)];
    ssdata          = [ssdata;sdata(:)];
end
At                  = sparse(ssi,ssj,ssdata,ndelta,m);
% fprintf('Done.\n')
end