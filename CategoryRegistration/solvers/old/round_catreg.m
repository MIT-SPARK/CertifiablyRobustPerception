function [R,t,c,theta] = round_catreg(X,N,K,opt)
if nargin < 4
    opt = 1;
end
X     = X{1};
[V,~] = sorteig(X);
%% take the opt-th eigenvector
nropt = length(opt);
x     = V(:,opt);
R     = zeros(3,3,nropt);
t     = zeros(3,nropt);
c     = zeros(K,nropt);
theta = zeros(N,nropt);
for k = 1:nropt
    xk    = x(:,k);
    xk    = xk/xk(1);
    %% round that eigenvector
    Rk    = reshape(xk(2:10),3,3);
    Rk    = project2SO3(Rk);
    tk    = xk(11:13);
    ck    = xk(14:13+K);
    thetak= sign( xk(14+10*K:13+10*K+N) );
    thetak(thetak==0) = 1;
    
    R(:,:,k)     = Rk;
    t(:,k)       = tk;
    c(:,k)       = ck;
    theta(:,k)   = thetak;
end
end