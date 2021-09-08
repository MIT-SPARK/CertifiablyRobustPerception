function S  = ape_dual_from_chordal_dual(cS)
%% Assemble a dense dual variable from the solution of a chordal relaxation
%% Heng Yang, July 01, 2021
N      = (length(cS)-1)/4;
n      = 13*N + 13;
%% moment blk
S      = zeros(n,n);
for i = 1:N
    cSi     = cS{4*i-2};
    idx     = [(1:13)';13+i;13+N+blkIndices(i,12)];
    S(idx,idx) = S(idx,idx) + cSi;
end
idx         = (1:13)';
S(idx,idx)  = S(idx,idx) + cS{1};

%% first PSD sub blk
S1     = zeros(1+N,1+N);
for i = 1:N
    cSi     = cS{4*i-1};
    idx     = [1;i+1];
    S1(idx,idx) = S1(idx,idx) + cSi;
end

%% second PSD sub blk
S2     = zeros(1+N,1+N);
for i = 1:N
    cSi     = cS{4*i};
    idx     = [1;i+1];
    S2(idx,idx) = S2(idx,idx) + cSi;
end

%% third PSD sub blk
S3     = zeros(1+N,1+N);
for i = 1:N
    cSi     = cS{4*i+1};
    idx     = [1;i+1];
    S3(idx,idx) = S3(idx,idx) + cSi;
end

S   = {S;S1;S2;S3};
end