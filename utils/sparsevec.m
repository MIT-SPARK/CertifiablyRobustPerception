function Avect   = sparsevec(blk,Acell)
%% Implement vec of a cell of sparse matrices
%% This is much faster than calling vec(Acell{i})
%% The generated Avect is typically used for CDCS

n           = blk{1,2};
m           = length(Acell);
nsq         = n^2;

ssi         = [];
ssj         = [];
ssdata      = [];
% fprintf('vec a cell of sparse matrices, progress ...')
for i = 1:m
    % if rem(i,10000) == 1
    %     fprintf('%d/%d ',i,m);
    % end
    Ai              = Acell{i};
    [ii,jj,val]     = find(Ai);
    
    si              = n * (jj-1) + ii;
    sj              = ones(length(si),1) * i;
    sdata           = val;
    
    ssi             = [ssi;si(:)];
    ssj             = [ssj;sj(:)];
    ssdata          = [ssdata;sdata(:)];
end
Avect               = sparse(ssi,ssj,ssdata,nsq,m);
% fprintf('Done.\n')
end