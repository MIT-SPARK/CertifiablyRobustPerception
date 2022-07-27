function [Xmat,pobjround,info] = local_search_quasar(Zmat,C,rrPar,opt,roundonly)
%% Round and refine for QUASAR SDP
%% Input:
%% Zmat: SDP primal iterate
%% C: cost matrix of the SDP
%% rrPar: a structure where you can pass any data necessary for local search
%% opt: indices of eigenvectors for rounding
%% roundonly: if true, then only round a solution without NLP
%% Output:
%% Xmat: a rank-one SDP iterate
%% pobjround: primal SDP cost attained by Xmat

if nargin < 5
    roundonly = false;
end
if nargin < 4
    % default round the first two eigenvectors
    opt = [1,2]; 
end

if roundonly
    [q,theta] = quasar_round(Zmat,[1]);
    xtld      = kron([1;theta],q);
    pobjround = xtld' * C{1} * xtld;
    Xmat      = {xtld * xtld'};
else
    [q,theta] = quasar_round(Zmat,opt);
    Zmatround = {};
    pobjround = zeros(length(opt),1);
    for i = 1:length(opt)
        [Ztmp,~,~,pobjtmp] = quasar_manopt(C,q(:,i),theta(:,i));
        Zmatround{end+1} = Ztmp;
        pobjround(i) = pobjtmp;
    end
    pobjs            = pobjround;
    [pobjround,idx]  = min(pobjround);
    if pobjround == inf
        fprintf('        NLP fails to find a good solution, return a rounded solution only.\n');
        qopt         = q(:,1); thetaopt = theta(:,1);
        xopt         = kron([1;thetaopt],qopt);
        Xmat         = {xopt * xopt'};    
        nlpsuccess   = false;
    else
        Xmat         = Zmatround{idx};
        nlpsuccess   = true;
    end
end

if nargout > 2
    info.minidx     = idx;
    info.nlpsuccess = nlpsuccess;
    info.pobjs      = pobjs;
    info.diffpobj   = pobjs(1) - pobjround;
end
end


function [q,theta] = quasar_round(X,opt)
if nargin < 2
    opt = 1;
end
[V,D] = eig(X{1});
[~,I] = sort(diag(D),'descend');
V     = V(:,I);
%% take the opt-th eigenvector
x     = V(:,opt);
q     = [];
theta = [];
for k = 1:size(x,2)
    xk    = x(:,k);
    %% round that eigenvector
    n     = length(xk);
    N     = n/4 - 1;
    qk    = xk(blkIndices(1,4));
    qk    = qk/norm(qk);
    thetak= zeros(N,1);
    for i = 1:N
        inprod = qk' * xk(blkIndices(i+1,4));
        if inprod > 0
            thetak(i) = 1;
        else
            thetak(i) = -1;
        end
    end
    q     = [q,qk];
    theta = [theta,thetak];
end
end

function [Xmat,qnew,thetanew,fopt] = quasar_manopt(C,q,theta)
N           = length(theta);
s           = spherefactory(4,1);
o           = obliquefactory(1,N);
elements.A  = s;
elements.B  = o;
manifold    = productmanifold(elements);

problem.M   = manifold;

warning('off', 'manopt:getHessian:approx') 
% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(x) quasar_cost(x,C);
problem.egrad = @(x) quasar_egrad(x,C);

% Numerically check gradient consistency (optional).
% checkgradient(problem);
% Solve.
x0.A                 = q;
x0.B                 = theta';
options.verbosity    = 0;
options.tolgradnorm  = 1e-6;
[xopt, fopt, output, options] = trustregions(problem,x0,options);
qopt = xopt.A;
thetaopt = xopt.B';
constraintviolationq          = abs(norm(qopt) - 1);
constraintviolationtheta      = max(abs(thetaopt.^2 - 1));
firstorderopt                 = output(end).gradnorm;
fprintf('        MANOPT: itr: %3d, constraint violation: %3.2e, %3.2e, gradnorm: %3.2e, cost: %3.8e.\n',...
    length(output),constraintviolationq,constraintviolationtheta,firstorderopt,fopt);
if max(constraintviolationq,constraintviolationtheta) < 1e-8 ...
        && firstorderopt < 1e-6
    % Do nothing
else
    fopt    = inf;
end
qnew = qopt/norm(qopt);
thetanew = sign(thetaopt);
vnew     = kron([1;thetanew],qnew);
Xmat     = {vnew * vnew'};
end


function f = quasar_cost(x,C)
q       = x.A;
theta   = x.B;
v       = kron([1;theta(:)],q);
f       = v'*C{1}*v;
end

function g = quasar_egrad(x,C)
q       = x.A;
theta   = x.B;
N       = length(theta);
v       = kron([1;theta(:)],q);

gv      = 2 * v' * C{1}; % 1 x n

vdq     = kron([1;theta(:)],speye(4)); % n x 4
gq      = (gv * vdq)'; % 4 x 1

vdtheta = [sparse(4,N);kron(speye(N),q)]; % n x N
gtheta  = (gv * vdtheta); % 1 x N

g.A     = gq;
g.B     = gtheta;
end
