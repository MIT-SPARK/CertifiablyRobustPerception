function SDP = chordal_relax_single_rotation_averaging(problem)
%% Apply a chordal sparse second-order relaxation to single rotation averaging
%% Generate SDP relaxations with multiple smaller blocks
%% This relaxation is likely to be loose than the single block relaxation
%% Depending on multivariate polynomial package SPOT
%% Heng Yang, June 29, 2021

fprintf('\n===================================================================')
fprintf('\nApplying Chordal SDP relaxation to single rotation averaging')
fprintf('\n===================================================================\n')
t0              = tic;

%% define POP variables
N               = problem.N;
noiseBoundSq    = problem.noiseBoundSq;
R_measurements  = problem.R_measurements;

nrPrimalVars    = 9 + N;
p               = msspoly('p',nrPrimalVars);
r               = p(1:9); 
col1            = r(1:3);
col2            = r(4:6);
col3            = r(7:9);
theta           = p(10:nrPrimalVars);
x               = r;

%% compute the cost function
residuals = [];
for i = 1:N 
    ri               = reshape(R_measurements(:,:,i),[9,1]);
    residuals        = [residuals; (6-2*ri'*r) / noiseBoundSq];
end
f_cost = [];
for i = 1:N 
    f_cost           = [f_cost; (1+theta(i))/2 * residuals(i) + (1-theta(i))/2 * 1.0];
end

%% Define the equality constraints
h_x = [1.0-col1'*col1;...
       1.0-col2'*col2;...
       1.0-col3'*col3;... % columns unit length
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

%% Formulate the chordal sparse second-order relaxation
%% the 0-th block [1;x] * [1;x]'
basis0          = [1;x];
n0              = length(basis0);
basis_x0        = 1;

pop0            = [mykron(basis_x0,h_x);...
                   mykron(basis0,basis0)];
[~,degmat,coef_all] = decomp(pop0);
coef_all            = coef_all';
dim_loc0        = length(basis_x0) * length(h_x);
n0delta         = triangle_number(n0);
nterms          = size(degmat,1);   
m_mom0          = n0delta - nterms;

assert(m_mom0==0,'The zero-th blk should have 0 moment constraints.')

coef_mom    = coef_all(:,dim_loc0+1:end);
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

%% the 1-N blocks [1;x;theta(i);theta(i)*x]' * [1;x;theta(i);theta(i)*x]
Aall            = {};
ball            = [];
Call            = {};
A0append        = {};
for blkidx = 1:N
    basis       = [1;x;theta(blkidx);theta(blkidx)*x];
    n           = length(basis);
    basis_x     = monomials(theta(blkidx),1:2);
    basis_theta = monomials(x,0:2);
    
    pop = [mykron(basis_x,h_x);...
           mykron(basis_theta,h_theta(blkidx));...
           mykron(basis,basis);
           f_cost(blkidx)];
    [~,degmat,coef_all] = decomp(pop);
    coef_all            = coef_all';
    
    % monomials_mom = mono(mykron(basis,basis));
    % fprintf(' %d ... %d ...\n',size(degmat,1),length(monomials_mom));
    % assert(size(degmat,1) == length(monomials_mom),'monomials not consistent');
    
    dim_loc             = length(basis_x) * length(h_x) + length(basis_theta) * length(h_theta(blkidx));
    [Acell,b,C]         = poly2SDP(n,degmat,coef_all,dim_loc,'verbose',0);
    Acell(1)  = []; Acell = Acell'; b(1) = [];
    
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
    Call             = [Call;{C}];
end

%% Convert to standard SDPT3 format
blk         = cell(N+1,2);
blk{1,1}    = 's'; blk{1,2} = n0;
for i = 1:N
    blk{i+1,1} = 's';
    blk{i+1,2} = n;
end
b           = [b0;ball];

A0t     = sparsesvec(blk(1,:),A0);
for i = 1:N
    A0t = [A0t,...
           sparse(n0delta,length(Aall{i})-length(A0append{i})),...
           sparsesvec(blk(1,:),A0append{i})];
end
ndelta  = triangle_number(n);
At    = {A0t};
for i = 1:N
    Ait = [sparse(ndelta,length(b0)),...
           sparse(ndelta,(i-1)*length(Aall{i})),...
           sparsesvec(blk(i+1,:),Aall{i}),...
           sparse(ndelta,(N-i)*length(Aall{i}))];
    At  = [At;{Ait}];
end

SDP.blk = blk;
SDP.At  = At;
SDP.n   = n;
SDP.m   = length(b);
SDP.C   = [{sparse(n0,n0)};Call];
SDP.b   = b;

tf    = toc(t0);
fprintf('\nDone in %g seconds.\n',tf);
fprintf('===================================================================\n')

%% Convert to Sedumi format
fprintf('Convert to Sedumi format ...\n')
t0    = tic;
sK.s  = [n0];
for i = 1:N 
    sK.s    = [sK.s,n];
end

A0t     = sparsevec(blk(1,:),A0);
n0sq    = n0^2;
for i = 1:N
    A0t = [A0t,...
           sparse(n0sq,length(Aall{i})-length(A0append{i})),...
           sparsevec(blk(1,:),A0append{i})];
end
nsq   = n^2;
At    = {A0t};
for i = 1:N
    Ait = [sparse(nsq,length(b0)),...
           sparse(nsq,(i-1)*length(Aall{i})),...
           sparsevec(blk(i+1,:),Aall{i}),...
           sparse(nsq,(N-i)*length(Aall{i}))];
    At  = [At;{Ait}];
end

sdata.K     = sK;
sdata.At    = cat(1,At{:});
sdata.b     = b;

sc          = [];
for i = 1:length(SDP.C)
    sc      = [sc;sparsevec(blk(i,:),SDP.C(i))];
end
sdata.c     = sc;

SDP.sedumi  = sdata;
tf    = toc(t0);
fprintf('Done in %g seconds.\n',tf);
fprintf('===================================================================\n')

end