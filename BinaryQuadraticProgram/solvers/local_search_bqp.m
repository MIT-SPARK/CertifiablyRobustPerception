function [Xr,pobj] = local_search_bqp(Xmat,Cmat,par,opt,roundonly)
if nargin < 5
    roundonly = false;
end

if nargin < 4
    opt = [1,2]; % default round the first two eigenvectors
end

d    = par.d;
v    = par.v;

if roundonly
    x       = bqp_round(Xmat,[1],d);
    v_est   = v.coefficient * (prod(x.^v.degmat,1))';
    Xr      = v_est*v_est';
    pobj    = v_est'*Cmat{1}*v_est;
else
    nhypo   = length(opt);
    x       = bqp_round(Xmat,opt,d);
    xopt    = zeros(d,nhypo);
    fopt    = zeros(nhypo,1);
    for i = 1:nhypo
        [xopti,fopti] = bqp_nlp(x(:,i),Cmat,par);
        xopt(:,i)     = xopti;
        fopt(i)       = fopti;
    end
    
    [pobj,idxmin]    = min(fopt);
    if pobj == inf
        fprintf('        NLP fails to find a good solution, return a rounded solution only.\n');
        xopt         = x(:,1);
        vopt         = v.coefficient * (prod(xopt.^v.degmat,1))';
        pobj         = vopt'*Cmat{1}*vopt;
        Xr           = {vopt*vopt'};
    else
        xopt         = xopt(:,idxmin);
        vopt         = v.coefficient * (prod(xopt.^v.degmat,1))';
        Xr           = {vopt*vopt'};
    end
end
end


function x = bqp_round(Xmat,opt,d)
nhypo       = length(opt);
Xmat        = Xmat{1};
X           = Xmat(1:d+1,1:d+1); % order 0 and 1 moment matrix
[V,D]       = eig(X);
[~,idx]     = sort(diag(D),'descend');
V           = V(:,idx);
xraw        = V(:,opt);
x           = zeros(d,nhypo);
for i = 1:nhypo
    xrawi   = xraw(:,i);
    xrawi   = xrawi/xrawi(1); % normalize first entry to be 1
    xi      = xrawi(2:end);
    
    xi      = sign(xi);
    xi(xi==0) = 1.0;
    x(:,i)  = xi;
end
end


function [xopt,fopt] = bqp_nlp(x0,Cmat,POP)
%% Nonlinear programming method for solving random POP
%% implemented using matlab fmincon
f           = POP.f; % cost function
J           = POP.J; % Jacobian
d           = length(x0);

cost        = @(x) randpop_cost(x,Cmat,f,J);
nonlincon   = @(x) binary_nonlincon(x);
options     = optimoptions('fmincon',...
               'Algorithm','interior-point',... % interior-point, sqp
               'Display','off',...
               'CheckGradients',false,...
               'MaxIterations',1000,...
               'SpecifyObjectiveGradient',true,...
               'SpecifyConstraintGradient',true);

[xopt,fopt,exitflag,output] = fmincon(cost,x0,[],[],[],[],-1*ones(d,1),ones(d,1),nonlincon,options);
%% Check if solution is actually good, if bad, return empty
fprintf('        NLP: exitflag: %d, constraint violation: %3.2e, first-order opt: %3.2e, cost: %3.8e.\n',exitflag,output.constrviolation,output.firstorderopt,fopt);
if exitflag > 0 || (output.constrviolation < 1e-8 && output.firstorderopt < 1e-3)
    % Do nothing
else
    fopt = inf; % this is not a good solution, return fopt = Inf;
end

end


function [f,g] = randpop_cost(x,Cmat,f_pop,J_pop)
f       = f_pop.coefficient * (prod(x.^f_pop.degmat,1))';
g       = J_pop.coefficient * (prod(x.^J_pop.degmat,1))';
end


function [c,ceq,gc,gceq] = binary_nonlincon(x)
c       = [];
gc      = [];
d       = length(x);
ceq     = x.^2 - 1;
gceq    = 2 * sparse(1:d,1:d,x,d,d);
end
