function S = catreg_dual_from_chordal_dual_v2(cS,N,K)
%% Assemble a dense dual variable from the solution of a chordal relaxation
%% Heng Yang, July 05, 2021
n      = (13+K)*(N+1);
nrIneq = 2+K;
%% moment blk
S      = zeros(n,n);
for i = 1:N
    idx = [(1:13+K)';...
           13+K+i;...
           13+K+N+blkIndices(i,12+K)];
    S(idx,idx) = S(idx,idx) + cS{(nrIneq+1)*i-nrIneq+1};
end
idx        = (1:13+K)';
S(idx,idx) = S(idx,idx) + cS{1};
S          = {S};

for k = 1:nrIneq
    Sk      = zeros(1+N,1+N);
    for i = 1:N
        idx = [1;i+1];
        Sk(idx,idx)  = Sk(idx,idx) + cS{(nrIneq+1)*i-nrIneq+k+1};
    end
    S = [S;{Sk}];
end
end