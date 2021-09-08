function [Xmat,pobjround,info] = local_search_pcr_v4(Zmat,C,rrPar,opt,roundonly)
%% local search for point cloud registration using MANOPT
%% used as a subroutine for STRIDE
%% Heng Yang
%% June 25, 2021

if nargin < 5
    roundonly = false;
end
if nargin < 4
    % default round the first two eigenvectors
    opt = [1,2]; 
end

tBound        = rrPar.translationBound;
blk           = rrPar.blk;

if roundonly
    [R,t,theta] = round_pcr_v4(Zmat,tBound,[1]);
    xtld        = lift_pcr_v4(R(:),t,theta,tBound);
    Xmat        = {xtld{1}*xtld{1}';xtld{2}*xtld{2}'};
    pobjround   = blktrace(blk,Xmat,C);
else
    [R,t,theta] = round_pcr_v4(Zmat,tBound,opt);
    Zmatround   = {};
    pobjround   = zeros(length(opt),1);
    for i = 1:length(opt)
        Ri              = squeeze(R(:,:,i));
        ti              = squeeze(t(:,i));
        thetai          = squeeze(theta(:,i));
        [Ztmp,pobjtmp]  = nlp_pcr(C,rrPar,Ri,ti,thetai);
        Zmatround{end+1}   = Ztmp;
        pobjround(i)       = pobjtmp;
    end
    pobjs            = pobjround;
    [pobjround,idx]  = min(pobjround);
    if pobjround == inf
        fprintf('        NLP fails to find a good solution, return a rounded solution only.\n');
        nlpsuccess   = false;
        Ropt         = squeeze(R(:,:,1));
        ropt         = Ropt(:);
        topt         = t(:,1);
        thetaopt     = theta(:,1);
        xopt         = lift_pcr_v4(ropt,topt,thetaopt,tBound);
        Xmat         = {xopt{1} * xopt{1}';xopt{2} * xopt{2}'};
        pobjround    = blktrace(blk,Xmat,C);          
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
    
    
function f = pcr_cost(x,C)
R       = x.A.R;
r       = R(:);
t       = x.A.t;
theta   = x.B;
theta   = theta(:);
v       = lift_pcr_v3(r,t,theta);
f       = v'*C*v;
end

function g = pcr_egrad(x,C)
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

function [Xmat,fopt] = nlp_pcr(Ccell,rrPar,R,t,theta)
blk         = rrPar.blk;
tBound      = rrPar.translationBound;
C           = Ccell{1}; % nxn matrix
N           = length(theta);
elements.A  = specialeuclideanfactory(3,1);
elements.B  = obliquefactory(1,N);
manifold    = productmanifold(elements);

problem.M   = manifold;

warning('off', 'manopt:getHessian:approx') 
% Define the problem cost function and its Euclidean gradient.
problem.cost  = @(x) pcr_cost(x,C);
problem.egrad = @(x) pcr_egrad(x,C);

% Numerically check gradient consistency (optional).
% checkgradient(problem);
% Solve.
x0.A.R               = R;
x0.A.t               = t;
x0.B                 = theta';
options.verbosity    = 0;
options.tolgradnorm  = 1e-6;
options.maxiter      = 100;
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
%     if norm(topt) > tBound
%         fopt = inf;
%         topt = topt/norm(topt) * tBound;
%     end
else
    fopt    = inf;
end
% Rnew = project2SO3(Ropt);
% rnew = Rnew(:);
% tnew = topt;
% thetanew = sign(thetaopt);
vnew     = lift_pcr_v4(Ropt(:),topt,thetaopt,tBound);
Xmat     = {vnew{1} * vnew{1}';vnew{2} * vnew{2}'};
end
        
            