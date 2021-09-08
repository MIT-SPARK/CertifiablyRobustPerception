function [Xmat,fopt,out] = nlp_catreg(Ccell,rrPar,R,t,c,theta)
blk         = rrPar.blk;
tBound      = rrPar.translationBound;
cBound      = rrPar.cBound;
N           = rrPar.N;
K           = rrPar.K;
C           = Ccell{1}; % nxn matrix
elements.A  = specialeuclideanfactory(3,1); % R and t
elements.B  = euclideanfactory(K,1); % shape params c
elements.C  = obliquefactory(1,N); % binary theta
manifold    = productmanifold(elements);

problem.M   = manifold;

warning('off', 'manopt:getHessian:approx') 
% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(x) catreg_cost(x,C,rrPar);
problem.egrad = @(x) catreg_egrad(x,C,rrPar);

% Numerically check gradient consistency (optional).
% checkgradient(problem);
% Solve.
x0.A.R               = R;
x0.A.t               = t;
x0.B                 = c;
x0.C                 = theta';
options.verbosity    = 0;
options.tolgradnorm  = 1e-6;
options.maxiter      = 1000;
% [xopt, fopt, output, options] = trustregions(problem,x0,options);
[xopt, fopt, output, options] = arc(problem,x0,options);
Ropt     = xopt.A.R;
topt     = xopt.A.t;
copt     = xopt.B;
thetaopt = xopt.C;
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

vnew     = lift_catreg(Ropt(:),topt,copt,thetaopt,cBound,tBound);
Xmat     = rank_one_lift(vnew);
if nargout > 2
    out.R    = Ropt;
    out.t    = topt;
    out.c    = copt;
    out.theta = thetaopt;
end
end


function f = catreg_cost(x,C,rrPar)
R       = x.A.R;
r       = R(:);
t       = x.A.t;
c       = x.B;
theta   = x.C;
theta   = theta(:);
v       = lift_catreg(r,t,c,theta,rrPar.translationBound,rrPar.cBound);
f       = v{1}'*C*v{1};
end

function g = catreg_egrad(x,C,rrPar)
R       = x.A.R;
r       = R(:);
t       = x.A.t;
c       = x.B;
theta   = x.C;
theta   = theta(:);
N       = rrPar.N;
K       = rrPar.K;
v       = lift_catreg(r,t,c,theta,rrPar.translationBound,rrPar.cBound);
v       = v{1};

gv      = 2 * v' * C; % 1 x n

vdr     =  [sparse(1,9);... % 1
            speye(9);... % r 
            sparse(3,9);... % t
            sparse(K,9);... % c
            kron(c,speye(9));... % c kron r
            sparse(N,9);... % theta
            kron(theta,[speye(9);sparse(3,9)]);... % theta kron [r;t]
            kron(theta,kron(c,speye(9)))]; % theta kron (c kron r)
gr      = (gv * vdr)'; % 9 x 1

vdt     = [ sparse(1,3);... % 1
            sparse(9,3);... % r
            speye(3);... % t
            sparse(K,3);... % c
            sparse(9*K,3);... % c kron r
            sparse(N,3);... % theta
            kron(theta,[sparse(9,3);speye(3)]);... % theta kron [r;t]
            sparse(9*K*N,3)]; % theta kron (c kron r)
gt      = (gv * vdt)'; % 3 x 1

vdc     = [sparse(1,K);... % 1
           sparse(9,K);... % r
           sparse(3,K);... % t
           speye(K);... % c
           kron(speye(K),r);... % c kron r
           sparse(N,K);... % theta
           sparse(12*N,K);... % theta kron [r;t]
           kron(theta,kron(speye(K),r))]; % theta kron (c kron r)
gc      = (gv * vdc)'; % K x 1

vdtheta = [sparse(13+10*K,N);... % 1,r,t,c,c kron r
            speye(N);... % theta
            kron(speye(N),[r;t]);... % theta kron [r;t]
            kron(speye(N),kron(c,r))]; % theta kron c kron r
gtheta  = (gv * vdtheta); % 1 x N

g.A.R   = reshape(gr,3,3);
g.A.t   = gt;
g.B     = gc;
g.C     = gtheta;
end