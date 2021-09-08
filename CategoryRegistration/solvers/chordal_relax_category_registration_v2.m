function SDP = chordal_relax_category_registration_v2(problem,varargin)
%% Apply a sparse, chordal, SECOND-order relaxation to category registration
%% Depending on multivariate polynomial package in SPOT
%% residual of v1: b(i) - R * (sum_k c_k ak(i)) - t
%% residual of v2: R*b(i) + t - sum_k c_k ak(i)
%% Heng Yang
%% July 06, 2021

params = inputParser;
params.CaseSensitive = false;

params.addParameter('checkMonomials',true, @(x) islogical(x));
params.addParameter('lambda',0.1, @(x) isscalar(x));

params.parse(varargin{:});

checkMonomials = params.Results.checkMonomials;
lambda         = params.Results.lambda;

fprintf('\n===================================================================')
fprintf('\nApplying Chordal SDP relaxation to category registration problem')
fprintf('\n===================================================================\n')
t0              = tic;

N               = problem.N;
K               = problem.K;
scene           = problem.scene;
shapes          = problem.shapes;
noiseBoundSq    = problem.noiseBoundSq;
tBound          = problem.translationBound;
tBoundSq        = tBound^2; % t'*t <= tBoundSq
cBoundSq        = problem.cBound^2; % should just be 1
barc2           = 1.0;

%% define POP variables
nrPrimalVars    = 9+3+K+N; % rotation: 9, translation: 3, binary: N, shape: K
p               = msspoly('p',nrPrimalVars);
r               = p(1:9);
R               = reshape(r,3,3); 
col1 = R(:,1); col2 = R(:,2); col3 = R(:,3);
t               = p(10:12);
c               = p(12+1:12+K);
theta           = p(12+K+1:nrPrimalVars);
x               = [r;t;c];

%% define cost function
shape           = combine_shapes(shapes,c);
residuals       = [];
for i = 1:N 
    distance            = R*scene(:,i) + t - shape(:,i);
    residuals           = [residuals; (distance' * distance) / noiseBoundSq];
end
f_cost = lambda * (c'*c); % regularization term
for i = 1:N
    f_cost = [f_cost;(1+theta(i))/2 * residuals(i) + (1-theta(i))/2 * barc2];
end

%% define constraints
h_r  = [1.0-col1'*col1;...
        1.0-col2'*col2;...
        1.0-col3'*col3;... % column unit length
        col1'*col2;...
        col2'*col3;...
        col3'*col1;... % colums orthogonal
        cross(col1,col2) - col3;...
        cross(col2,col3) - col1;...
        cross(col3,col1) - col2]; % columns righthandedness

h_theta = [];
for i = 1:N 
    h_theta =[h_theta; 1-theta(i)^2];
end

g_t = tBoundSq - t'*t; % Translation bounded
g_c = [cBoundSq - c'*c;c]; % nonnegative and bounded shape parameters
g_x = [g_t;g_c];

%% Formulate the chordal sparse second-order relaxation
%% the 0-th block [1;x] * [1;x]'
basis0          = [1;x];
n0              = length(basis0);
basis_x0        = get_multiplier_basis(x,basis0,h_r(1));
pop0            = [mykron(basis_x0,h_r);...
                   mykron(basis0,basis0);...
                   f_cost(1)];
[~,degmat,coef_all] = decomp(pop0);
coef_all            = coef_all';
dim_loc0        = length(basis_x0) * length(h_r);
n0delta         = triangle_number(n0);
nterms          = size(degmat,1);   
m_mom0          = n0delta - nterms;

assert(m_mom0==0,'The zero-th blk should have 0 moment constraints.')

coef_mom    = coef_all(:,dim_loc0+1:dim_loc0+n0^2);
coef_mom    = coef_mom';
B           = {};
B_normalize = {};

for i = 1:nterms
    [row,~,~]   = find(coef_mom(:,i));
    SDP_coli    = floor((row-1)./n0) + 1;
    SDP_rowi    = mod(row-1,n0) + 1;
    nnz         = length(SDP_rowi);
    
    Bi          = sparse(SDP_rowi,SDP_coli,ones(nnz,1),n0,n0);
    B{end+1}    = Bi;
    B_normalize{end+1} = Bi/nnz;
end

coef_loc0       = coef_all(:,1:dim_loc0);
A0_local        = {};

for i = 1:dim_loc0
    [rowi,~,vi] = find(coef_loc0(:,i));
    Ai      = sparse(n0,n0);
    for j   = 1:length(rowi)
        Ai  = Ai + vi(j) * B_normalize{rowi(j)};
    end
    A0_local = [A0_local;{Ai}];
end
A0_0    = sparse([1],[1],[1],n0,n0);
% The first block satisfies A0(X0) = b0;
A0      = [{A0_0};A0_local];
b0      = sparse(1,1,1,length(A0),1);

% Now build the cost matrix
coef_cost   = coef_all(:,dim_loc0+n0^2+1);
[row,~,v]   = find(coef_cost);
C           = sparse(n0,n0);
for i = 1:length(row)
    C       = C + v(i) * B_normalize{row(i)};
end

%% the 1-N blocks [1;x;theta(i);theta(i)*x] * [1;x;theta(i);theta(i)*x]'
%% Since there are (2+K) inequality constraints, it will generate (K+3)*N blocks
nrineq          = length(g_x);
Aall            = {};
Asuball         = {};
ball            = [];
Call            = {C};
A0append        = {};
for blkidx = 1:N
    basis       = [1;x;theta(blkidx);theta(blkidx)*x];
    n           = length(basis);
    basis_x     = get_multiplier_basis([x;theta(blkidx)],basis,h_r(1),0);
    basis_theta = get_multiplier_basis([x;theta(blkidx)],basis,h_theta(blkidx),0);
    basis_g     = [1;theta(blkidx)];
    
    out         = gen_chordal_subblk_catreg_v2(...
        basis,basis_x,h_r,basis_theta,h_theta(blkidx),g_x,basis_g,f_cost(1+blkidx));
    
    Acell       = out.A;
    Asub        = out.Asub;
    b           = out.b;
    C           = out.C;
    
    % add constraint that the top-left [1;x]*[1;x]' block is the same as
    % the 0-th block
    A0blk = {};
    for i = 1:n0
        for j = i:n0
            if i == j
                A0i  = sparse(i,j,-1,n0,n0);
                Ai   = sparse(i,j,1,n,n);
            else
                A0i  = sparse([i,j],[j,i],[-0.5,-0.5],n0,n0);
                Ai   = sparse([i,j],[j,i],[0.5,0.5],n,n);
            end
            A0blk    = [A0blk;{A0i}];
            Acell    = [Acell;{Ai}];
        end
    end
    
    ball             = [ball;b;sparse(n0delta,1)];
    A0append{end+1}  = A0blk;
    Aall{end+1}      = Acell;
    Call             = [Call;C];
    Asuball{end+1}   = Asub;
end

%% Convert to standard SDPT3 format
b           = [b0;ball];
blk         = cell( (nrineq+1)*N+1,2);
blk{1,1}    = 's'; blk{1,2} = n0;
n1          = out.blk{2,2};
n1delta     = triangle_number(n1);
step        = nrineq + 1;
for i = 1:N
    blk{step*i-nrineq+1,1} = 's';
    blk{step*i-nrineq+1,2} = n;
    for j = 1:nrineq
        blk{step*i-nrineq+j+1,1} = 's';
        blk{step*i-nrineq+j+1,2} = n1;
    end
end

A0t     = sparsesvec(blk(1,:),A0);
for i = 1:N
    A0t = [A0t,...
           sparse(n0delta,out.m),...
           sparsesvec(blk(1,:),A0append{i})];
end

ndelta  = triangle_number(n);
At      = {A0t};
for i = 1:N
    Ait = [sparse(ndelta,length(b0)),...
           sparse(ndelta,(i-1)*length(Aall{i})),...
           sparsesvec(blk(step*i-nrineq+1,:),Aall{i}),... % the moment matrix block
           sparse(ndelta,(N-i)*length(Aall{i}))];
    At  = [At;{Ait}];
    for j = 1:nrineq
        Aijt    = [sparse(n1delta,length(b0)),...
                    sparse(n1delta,(i-1)*length(Aall{i})),...
                    sparse(n1delta,out.m_mom+out.m_loc),... % the moment constraint
                    sparse(n1delta,(j-1)*n1delta),...
                    sparsesvec(blk(step*i-nrineq+j+1,:),Asuball{i}{j}),...
                    sparse(n1delta,(nrineq-j)*n1delta),...
                    sparse(n1delta,n0delta),...
                    sparse(n1delta,(N-i)*length(Aall{i}))];

        At  = [At;{Aijt}];
    end
end

SDP.blk = blk;
SDP.At  = At;
SDP.m   = length(b);
SDP.C   = Call;
SDP.b   = b;

tf    = toc(t0);
fprintf('\nDone in %g seconds.\n',tf);
fprintf('===================================================================\n')

%% Convert to Sedumi format
fprintf('Convert to Sedumi format ...\n')
t0    = tic;
sK.s  = [n0];
for i = 1:N
    sK.s        = [sK.s,n,n1*ones(1,nrineq)];
end

A0t     = sparsevec(blk(1,:),A0);
n0sq    = n0^2;
for i = 1:N
    A0t = [A0t,...
           sparse(n0sq,out.m),...
           sparsevec(blk(1,:),A0append{i})];
end

nsq     = n^2;
n1sq    = n1^2;
At      = {A0t};
for i = 1:N
    Ait = [sparse(nsq,length(b0)),...
           sparse(nsq,(i-1)*length(Aall{i})),...
           sparsevec(blk(step*i-nrineq+1,:),Aall{i}),... % the moment matrix block
           sparse(nsq,(N-i)*length(Aall{i}))];
    At  = [At;{Ait}];
    for j = 1:nrineq
        Aijt    = [sparse(n1sq,length(b0)),...
                    sparse(n1sq,(i-1)*length(Aall{i})),...
                    sparse(n1sq,out.m_mom+out.m_loc),... % the moment constraint
                    sparse(n1sq,(j-1)*n1delta),...
                    sparsevec(blk(step*i-nrineq+j+1,:),Asuball{i}{j}),...
                    sparse(n1sq,(nrineq-j)*n1delta),...
                    sparse(n1sq,n0delta),...
                    sparse(n1sq,(N-i)*length(Aall{i}))];

        At  = [At;{Aijt}];
    end
end

sdata.K     = sK;
sdata.At    = cat(1,At{:});
sdata.b     = b;

sc          = [];
for i = 1:length(Call)
    sc      = [sc;sparsevec(blk(i,:),Call(i))];
end
sdata.c     = sc;

SDP.sedumi   = sdata;


tf    = toc(t0);
fprintf('Done in %g seconds.\n',tf);
fprintf('===================================================================\n')


end