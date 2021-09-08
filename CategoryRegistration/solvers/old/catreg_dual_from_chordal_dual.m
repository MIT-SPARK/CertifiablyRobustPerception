function S = catreg_dual_from_chordal_dual(cS,N,K)
%% Assemble a dense dual variable from the solution of a chordal relaxation
%% Heng Yang, July 05, 2021
n      = 13+10*K+13*N+9*K*N;
%% moment blk
S      = zeros(n,n);
for i = 1:N 
    idx = [(1:13+10*K)';...
          13+10*K+i;...
          13+10*K+N+blkIndices(i,12);...
          13*N+10*K+13+blkIndices(i,9*K)];
    S(idx,idx) = S(idx,idx) + cS{1+K+2*i};
end
idx = (1:13+10*K)';
S(idx,idx) = S(idx,idx) + cS{1};

%% translation blk
St     = zeros(1+N,1+N);
for i = 1:N 
    idx = [1;i+1];
    St(idx,idx) = St(idx,idx) + cS{1+K+2*i+1};
end

%% shape params blk
Sc     = cS(2:K+2);

S      = [{S};{St};Sc];
end