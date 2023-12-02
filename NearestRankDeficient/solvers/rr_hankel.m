function [Xr,pobj] = rr_hankel(Xmat,Cmat,par,opt,roundonly)
%% Implement rounding and refine for nearest rank-deficient Hankel approximation
if nargin < 5
    roundonly = false;
end

if nargin < 4
    opt = [1,2]; % default round the first two eigenvectors
end

m          = par.m;
k          = par.k;
Xmat       = Xmat{1};
Cmat       = Cmat{1};

if roundonly
    [z,u]  = hankel_round(Xmat,[1],par);
    xtld   = kron([u;1],z);
    Xr     = xtld * xtld';
    pobj   = xtld' * Cmat * xtld;
else
    numhypo  = length(opt); % number of hypotheses
    [z,u]    = hankel_round(Xmat,opt,par);
    zopt     = zeros(m,numhypo);
    uopt     = zeros(k,numhypo);
    fopt     = zeros(numhypo,1);
    
    for i = 1:numhypo
        [zopti,uopti,fopti] = hankel_nlp(z(:,i),u(:,i),par);
        zopt(:,i)           = zopti;
        uopt(:,i)           = uopti;
        fopt(i)             = fopti;
    end
    
    [pobj,idxmin]           = min(fopt);
    
    if pobj == inf
        fprintf('        NLP fails to find a good solution, return a rounded solution only.\n');
        zr                      = z(:,1);
        ur                      = u(:,1);
        xtld                    = kron([ur;1],zr);
        Xr                      = {xtld * xtld'};
    else
        zr                      = zopt(:,idxmin);
        ur                      = uopt(:,idxmin);
        xtld                    = kron([ur;1],zr);
        Xr                      = {xtld * xtld'};
    end
    % fprintf('pobj by nlp: %3.2e, pobj by sdp: %3.2e.\n',pobj,xtld'*Cmat*xtld);
end

end


function [z,u] = hankel_round(Xmat,opt,par)
%% Rounding from moment matrix Xmat
m       = par.m;
k       = par.k;

[V,D]   = eig(Xmat);
[~,idx] = sort(diag(D),'descend');
V       = V(:,idx);

nhypo   = length(opt); % number of hypotheses
xraw    = V(:,opt); % (k+1)m by nhypo

z       = zeros(m,nhypo);
u       = zeros(k,nhypo);

for i = 1:nhypo
    xi  = xraw(:,i); % (k+1)m by 1
    zi  = xi(blkIndices(k+1,m)); % trailing m by 1 block is z
    ui  = zeros(k,1);
    for j = 1:k
        ui(j) = zi'*xi(blkIndices(j,m));
    end
    
    z(:,i) = zi/norm(zi); % make sure z is unit norm
    u(:,i) = ui;
end
end


function [zopt,uopt,fopt] = hankel_nlp(z0,u0,par)
%% Solve z and u from nonlinear programming with initial guess
theta   = par.theta;
S       = par.S;
m       = par.m;
k       = par.k;

cost       = @(x) hankel_cost(x,theta,m);
nonlincon  = @(x) hankel_nonlincon(x,S,m,k);

x0      = [z0;u0];

options = optimoptions('fmincon',...
                       'Algorithm','interior-point',...
                       'Display','off',...
                       'CheckGradients',false,...
                       'MaxIterations',1000,...
                       'SpecifyObjectiveGradient',true,...
                       'SpecifyConstraintGradient',true,...
                       'OptimalityTolerance',1e-12,...
                       'StepTolerance',1e-32,...
                       'ConstraintTolerance',1e-12);

[xopt,fopt,exitflag,output] = fmincon(cost,x0,[],[],[],[],[],[],nonlincon,options);

%% Check if solution is actually good, if bad, return empty
fprintf('        NLP: exitflag: %d, itr: %4d, constraint violation: %3.2e, first-order opt: %3.2e, cost: %3.8e.\n',...
    exitflag,output.iterations,output.constrviolation,output.firstorderopt,fopt);

if exitflag > 0 || (output.constrviolation < 1e-8 && output.firstorderopt < 1e-3)
    % Do nothing
else
    fopt = inf; % this is not a good solution, return fopt = Inf;
end    

zopt        = xopt(1:m);
uopt        = xopt(m+1:end);
end


function [f,g] = hankel_cost(x,theta,m)
% x = (z || u)
% z       = x(1:m);
u       = x(m+1:end);

guhalf  = u - theta;
f       = norm(guhalf)^2;

g       = [zeros(m,1);2*guhalf]; % (m+k) by 1
end

function [c,ceq,gc,gceq] = hankel_nonlincon(x,S,m,k)
% S = (S1 || S2 || ... || Sk || S0) \in (k+1)m x n
St      = S';
St1     = St(:,1:k*m); % n by km
% St1 = (S1', S2', ... Sk')

c       = [];
gc      = [];

z       = x(1:m);
u       = x(m+1:end);

Im      = speye(m);
utld    = [u;1];
Su      = kron(utld',Im) * S; % m by n
Sut     = Su'; % n by m

ceq     = [z'*z - 1; Sut*z]; % (n+1) by 1

Ik      = speye(k);
gceq    = [2*z',sparse(1,k);
           Sut,St1*kron(Ik,z)]; % (n+1) by (m+k)
gceq    = gceq';
end