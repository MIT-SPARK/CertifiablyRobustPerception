function [Xmat,pobjround,info] = local_search_sra(Zmat,C,rrPar,opt,roundonly)
%% local search for single rotation averaging 
%% used as a subroutine for STRIDE
if nargin < 5
    roundonly = false;
end
if nargin < 4
    % default round the first two eigenvectors
    opt = [1,2]; 
end

if roundonly
    [R,theta] = round_sra(Zmat,[1]);
    xtld      = lift_sra(R(:),theta);
    pobjround = xtld{1}' * C{1} * xtld{1};
    Xmat      = rank_one_lift(xtld);
else
    [R,theta] = round_sra(Zmat,opt);
    Zmatround = {};
    pobjround = zeros(length(opt),1);
    for i = 1:length(opt)
        Ri              = squeeze(R(:,:,i));
        thetai          = theta(:,i);
        [Ztmp,~,~,pobjtmp] = nlp_sra(C,Ri,thetai);
        Zmatround{end+1} = Ztmp;
        pobjround(i) = pobjtmp;
    end
    pobjs            = pobjround;
    [pobjround,idx]  = min(pobjround);
    if pobjround == inf
        fprintf('        NLP fails to find a good solution, return a rounded solution only.\n');
        Ropt         = squeeze(R(:,:,1));
        ropt         = Ropt(:);
        thetaopt     = theta(:,1);
        xopt         = lift_sra(ropt,thetaopt);
        Xmat         = rank_one_lift(xopt);
        pobjround    = xopt{1}'*C{1}*xopt{1};
        nlpsuccess   = false;
    else
        Xmat = Zmatround{idx};
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


function f = sra_cost(x,C)
R       = x.A;
r       = R(:);
theta   = x.B;
theta   = theta(:);
v       = lift_sra(r,theta);
f       = v{1}'*C*v{1};
end

function g = sra_egrad(x,C)
R       = x.A;
r       = R(:);
theta   = x.B;
theta   = theta(:);
N       = length(theta);
v       = lift_sra(r,theta);
v       = v{1};

gv      = 2 * v' * C; % 1 x n

vdr     = [sparse(1,9);...
           speye(9);...
           sparse(N,9);...
           kron(theta,speye(9))]; % n x 9
gr      = (gv * vdr)'; % 9 x 1

vdtheta = [sparse(10,N);...
           speye(N);...
           kron(speye(N),r)]; % n x N
gtheta  = (gv * vdtheta); % 1 x N

g.A     = reshape(gr,3,3);
g.B     = gtheta;
end

function [Xmat,Ropt,thetaopt,fopt] = nlp_sra(C,R,theta)
N           = length(theta);
elements.A  = rotationsfactory(3,1);
elements.B  = obliquefactory(1,N);
manifold    = productmanifold(elements);

problem.M   = manifold;

warning('off', 'manopt:getHessian:approx') 
% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(x) sra_cost(x,C{1});
problem.egrad = @(x) sra_egrad(x,C{1});

% Numerically check gradient consistency (optional).
% checkgradient(problem);
% Solve.
x0.A                 = R;
x0.B                 = theta';
options.verbosity    = 0;
options.tolgradnorm  = 1e-6;
[xopt, fopt, output, options] = trustregions(problem,x0,options);
Ropt = xopt.A;
thetaopt = xopt.B';
constraintviolationR          = norm(Ropt*Ropt'-eye(3),'fro');
constraintviolationtheta      = max(abs(thetaopt.^2 - 1));
firstorderopt                 = output(end).gradnorm;
fprintf('        MANOPT: itr: %3d, constraint violation: %3.2e, %3.2e, gradnorm: %3.2e, cost: %3.8e.\n',...
    length(output),constraintviolationR,constraintviolationtheta,firstorderopt,fopt);
if max(constraintviolationR,constraintviolationtheta) < 1e-8 ...
        && firstorderopt < 1e-6
    % Do nothing
else
    fopt    = inf;
end
vnew     = lift_sra(Ropt(:),thetaopt(:));
Xmat     = rank_one_lift(vnew);
end

