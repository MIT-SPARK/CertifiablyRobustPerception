function S = sra_dual_from_chordal_dual(cS)
%% Assemble a dense dual variable from the solution of a chordal relaxation
%% Heng Yang, June 29, 2021
N      = length(cS)-1;
n      = 10*N + 10;
S      = zeros(n,n);
for i = 1:N
    cSi     = cS{i+1};
    idx     = [(1:10)';10+i;10+N+blkIndices(i,9)];
    S(idx,idx) = S(idx,idx) + cSi;
end
idx         = (1:10)';
S(idx,idx)  = S(idx,idx) + cS{1};
S      = {S};
end