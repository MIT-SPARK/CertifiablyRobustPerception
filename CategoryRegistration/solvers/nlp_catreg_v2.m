function [Xmat,fopt,out] = nlp_catreg_v2(Ccell,rrPar,R,t,c,theta)
blk         = rrPar.blk;
tBound      = rrPar.translationBound;
cBound      = rrPar.cBound;
N           = rrPar.N;
K           = rrPar.K;
C           = Ccell{1}; % nxn matrix
elements.A  = specialeuclideanfactory(3,1); % R and t
% elements.B  = euclideanfactory(K,1); % shape params c positivefactory(m, n)
elements.B  = positivefactory(K,1);
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
options.maxiter      = 100;
[xopt, fopt, output, options] = trustregions(problem,x0,options);
% [xopt, fopt, output, options] = arc(problem,x0,options);
Ropt     = xopt.A.R;
topt     = xopt.A.t;
copt     = xopt.B;
thetaopt = xopt.C;
constraintviolationR          = norm(Ropt*Ropt'-eye(3),'fro');
constraintviolationtheta      = max(abs(thetaopt.^2 - 1));
constraintviolationc          = max(max(-copt,0));
firstorderopt                 = output(end).gradnorm;
fprintf('        MANOPT: itr: %3d, constraint violation: %3.2e, %3.2e, %3.2e, gradnorm: %3.2e, cost: %3.8e.\n',...
    length(output),constraintviolationR,constraintviolationtheta,constraintviolationc,firstorderopt,fopt);
if max([constraintviolationR,constraintviolationtheta,constraintviolationc]) < 1e-6 ...
        && firstorderopt < 1e-1
    % Do nothing
    % make sure c_est is norm bounded
    if norm(copt) >= cBound
        copt = copt / norm(copt) * cBound;
    end
    % make sure t_est is norm bounded
    if norm(topt) >= tBound
        topt = topt / norm(topt) * tBound;
    end
    xnew.A.R  = Ropt;
    xnew.A.t  = topt;
    xnew.B    = copt;
    xnew.C    = thetaopt;
    
    fopt      = catreg_cost(xnew,C,rrPar);
    
else
    fopt    = inf;
end

vnew     = lift_catreg_v2(Ropt(:),topt,copt,thetaopt,cBound,tBound);
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
v       = lift_catreg_v2(r,t,c,theta,rrPar.translationBound,rrPar.cBound);
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
v       = lift_catreg_v2(r,t,c,theta,rrPar.translationBound,rrPar.cBound);
v       = v{1};

gv      = 2 * v' * C; % 1 x n

vdr     =  [sparse(1,9);... % 1
            speye(9);... % r 
            sparse(3,9);... % t
            sparse(K,9);... % c
            sparse(N,9);... % theta
            kron(theta,[speye(9);sparse(3,9);sparse(K,9)])]; % theta kron [r;t;c]
gr      = (gv * vdr)'; % 9 x 1

vdt     = [ sparse(1,3);... % 1
            sparse(9,3);... % r
            speye(3);... % t
            sparse(K,3);... % c
            sparse(N,3);... % theta
            kron(theta,[sparse(9,3);speye(3);sparse(K,3)])]; % theta kron [r;t;c]

gt      = (gv * vdt)'; % 3 x 1

vdc     = [sparse(1,K);... % 1
            sparse(9,K);... % r
            sparse(3,K);... % t
            speye(K);... % c
            sparse(N,K);... % theta
            kron(theta,[sparse(9,K);sparse(3,K);speye(K)])]; % theta kron [r;t;c]
gc      = (gv * vdc)'; % K x 1

vdtheta = [sparse(13+K,N);... % 1,r,t,c
            speye(N);... % theta
            kron(speye(N),[r;t;c])]; % theta kron [r;t]
gtheta  = (gv * vdtheta); % 1 x N

g.A.R   = reshape(gr,3,3);
g.A.t   = gt;
g.B     = gc;
g.C     = gtheta;
end