function S  = pcr_dual_from_chordal_dual(cS)
%% Assemble a dense dual variable from the solution of a chordal relaxation
%% Heng Yang, June 29, 2021
N      = (length(cS)-1)/2;
n      = 13*N + 13;
S      = zeros(n,n);
for i = 1:N
    cSi     = cS{2*i};
    idx     = [(1:13)';13+i;13+N+blkIndices(i,12)];
    S(idx,idx) = S(idx,idx) + cSi;
end
idx         = (1:13)';
S(idx,idx)  = S(idx,idx) + cS{1};

S1     = zeros(1+N,1+N);
for i = 1:N
    cSi     = cS{2*i+1};
    idx     = [1;i+1];
    S1(idx,idx) = S1(idx,idx) + cSi;
end

S   = {S;S1};
end