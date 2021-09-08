function [R,theta] = round_sra(X,opt)
%% Given the moment matrix X, round R and theta from X
if nargin < 2
    opt = 1;
end
X     = X{1};
n     = size(X,1);
N     = round((n-10)/10);
[V,D] = eig(X);
[~,I] = sort(diag(D),'descend');
V     = V(:,I);
%% take the opt-th eigenvector
nropt = length(opt);
x     = V(:,opt);
R     = zeros(3,3,nropt);
theta = zeros(N,nropt);
for k = 1:nropt
    xk    = x(:,k);
    xk    = xk/xk(1);
    %% round that eigenvector
    Rk    = reshape( xk(1+blkIndices(1,9)),3,3 );
    Rk    = project2SO3(Rk);
    
    thetak= sign( xk(10+1:10+N) );
    thetak(thetak==0) = 1;
    
    R(:,:,k)     = Rk;
    theta(:,k)   = thetak;
end
end