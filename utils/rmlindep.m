function [SDP,idxkpt] = rmlindep(SDP,tol)
%% remove linearly dependent constraints
if nargin < 2
    tol = 1e-5;
end

m       = length(SDP.b);
AAt     = sparse(m,m);
blk     = SDP.blk;
blknum  = size(blk,1);
At      = SDP.At;

for i = 1:blknum
    AAt = AAt + At{i}'*At{i};
end

AAt     = AAt + 1e-12 * speye(m,m);
R       = chol(AAt);
dR      = diag(R);
idxkp   = (abs(dR) > tol);

m       = full(sum(idxkp));

% fprintf('\n\nKeep %d out of %d constraints.\n',m,length(SDP.b));

Atkp    = {};
for i = 1:blknum
    Ati     = At{i};
    Ati     = Ati(:,idxkp);
    Atkp    = [Atkp;{Ati}];
end

SDP.b   = SDP.b(idxkp);
SDP.At  = Atkp;
SDP.m   = m;

if nargout > 1
    idxkpt = idxkp;
end

end