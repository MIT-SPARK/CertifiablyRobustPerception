function [res,fopt] = nlp_ape(Ccell,rrPar,R,t,theta)
blk         = rrPar.blk;
tBound      = rrPar.translationBound;
dBound      = rrPar.depthBound;
FOV         = rrPar.FOV;
C           = Ccell{1};
N           = length(theta);
elements.A  = specialeuclideanfactory(3,1);
elements.B  = obliquefactory(1,N);
manifold    = productmanifold(elements);

problem.M   = manifold;

warning('off', 'manopt:getHessian:approx') 
% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(x) ape_cost(x,C);
problem.egrad = @(x) ape_egrad(x,C);

% Numerically check gradient consistency (optional).
% checkgradient(problem);
% Solve.
x0.A.R               = R;
x0.A.t               = t;
x0.B                 = theta';
options.verbosity    = 0;
options.tolgradnorm  = 1e-6;
options.maxiter      = 1000;
[xopt, fopt, output, options] = trustregions(problem,x0,options);
% [xopt, fopt, output, options] = arc(problem,x0,options);
Ropt     = xopt.A.R;
topt     = xopt.A.t;
thetaopt = xopt.B;
constraintviolationR          = norm(Ropt*Ropt'-eye(3),'fro');
constraintviolationtheta      = max(abs(thetaopt.^2 - 1));
firstorderopt                 = output(end).gradnorm;
fprintf('        MANOPT: itr: %3d, constraint violation: %3.2e, %3.2e, gradnorm: %3.2e, cost: %3.8e.\n',...
    length(output),constraintviolationR,constraintviolationtheta,firstorderopt,fopt);
if max([constraintviolationR,constraintviolationtheta]) < 1e-8 ...
        && firstorderopt < 1e-3
    % Do nothing
else
    fopt    = inf;
end
if ~check_translation(topt,tBound,FOV)
    fopt    = inf;
end
% Rnew = project2SO3(Ropt);
% rnew = Rnew(:);
% tnew = topt;
% thetanew = sign(thetaopt);
vnew     = lift_ape_v1(Ropt(:),topt,thetaopt,tBound,dBound,FOV);
Xmat     = {vnew{1} * vnew{1}';vnew{2} * vnew{2}';vnew{3} * vnew{3}';vnew{4} * vnew{4}'};
% fopt     = blktrace(blk,Xmat,Ccell);

res.X     = Xmat;
res.R     = Ropt;
res.t     = topt;
res.theta = thetaopt;
end

function f = ape_cost(x,C)
R       = x.A.R;
r       = R(:);
t       = x.A.t;
theta   = x.B;
theta   = theta(:);
v       = lift_pcr_v3(r,t,theta);
f       = v'*C*v;
end

function g = ape_egrad(x,C)
R       = x.A.R;
r       = R(:);
t       = x.A.t;
theta   = x.B;
theta   = theta(:);
N       = length(theta);
v       = lift_pcr_v3(r,t,theta);

gv      = 2 * v' * C; % 1 x n

vdr     = [sparse(1,9);...
            speye(9);...
            sparse(3,9);...
            sparse(N,9);...
            kron(theta,[speye(9);sparse(3,9)])]; % n x 9
gr      = (gv * vdr)'; % 9 x 1

vdt     = [ sparse(1,3);...
            sparse(9,3);...
            speye(3);...
            sparse(N,3);...
            kron(theta,[sparse(9,3);speye(3)])]; % n x 4
gt      = (gv * vdt)'; % 3 x 1

vdtheta = [sparse(13,N);...
            speye(N);...
            kron(speye(N),[r;t])]; % n x N
gtheta  = (gv * vdtheta); % 1 x N

g.A.R   = reshape(gr,3,3);
g.A.t   = gt;
g.B     = gtheta;
end