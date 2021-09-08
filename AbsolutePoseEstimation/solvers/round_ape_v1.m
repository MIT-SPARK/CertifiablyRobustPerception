function [R,t,theta] = round_ape_v1(Xcell,tBound,depthBound,opt)
    if nargin < 4
        opt = 1;
    end
    X     = Xcell{1};
    n     = size(X,1);
if mod(n-13,13) == 0
    N     = n/13 - 1;
elseif mod(n-40,13) == 0
    N     = (n-40)/13;
end
    [V,D] = eig(X);
    [~,I] = sort(diag(D),'descend');
    V     = V(:,I);
    %% take the opt-th eigenvector
    nropt = length(opt);
    x     = V(:,opt);
    R     = zeros(3,3,nropt);
    t     = zeros(3,nropt);
    theta = zeros(N,nropt);
    for k = 1:nropt
        xk    = x(:,k);
        xk    = xk/xk(1);
        %% round that eigenvector
        Rk    = reshape(xk(2:10),3,3);
        Rk    = project2SO3(Rk);
        tk    = xk(11:13);
        thetak= sign( xk(13+1:13+N) );
        thetak(thetak==0) = 1;
        
        R(:,:,k)     = Rk;
        t(:,k)       = tk;
        theta(:,k)   = thetak;
    end
end