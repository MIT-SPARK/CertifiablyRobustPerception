function SDP = chordal_relax_category_registration(problem,varargin)
%% Apply a sparse, chordal, third-order relaxation to category registration
%% Depending on multivariate polynomial package in SPOT
%% Bounded translation is modelled as an inequality constraint, which leads
%% to an extra PSD block in the semidefinite relaxation
%% Heng Yang
%% July 05, 2021

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
x               = [r;t];
cr              = mykron(c,r);

%% define cost function
shape           = combine_shapes(shapes,c);
residuals       = [];
for i = 1:N 
    distance            = scene(:,i) - R * shape(:,i) - t;
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

%% Formulate the chordal relaxation
%% First, v0 = [1;r;t;c;cr]
basis0          = [1;x;c;cr];
n0              = length(basis0);
basis_r0        = get_multiplier_basis(p(1:12+K),basis0,h_r(1));
basis_g_c       = [1;r];
n01             = length(basis_g_c);
basis_g_c_f     = mykron(basis_g_c,basis_g_c);

pop0            = [mykron(basis_r0,h_r);...
                   mykron(basis0,basis0);...
                   f_cost(1);...
                   mykron(g_c,basis_g_c_f)];
[~,degmat,coef_all] = decomp(pop0);
coef_all            = coef_all';
if checkMonomials
    fprintf('Checking consistency of monomials ...')
    time_check0   = tic;
    monomials_mom = mono(mykron(basis0,basis0));
    fprintf(' %d ... %d ...',size(degmat,1),length(monomials_mom));
    assert(size(degmat,1) == length(monomials_mom),'monomials not consistent');
    time_check    = toc(time_check0);
    fprintf('Done in %g seconds.\n',time_check);
end

n0delta         = triangle_number(n0);
n01delta        = triangle_number(n01);
dim_loc_eq      = length(basis_r0) * length(h_r);
dim_loc_ineq    = length(g_c) * n01delta;

%% generate standard SDP data from degmat and coefficients
nterms      = size(degmat,1);
m_mom       = n0delta - nterms;
m_loc       = dim_loc_eq;
m_loc_ineq  = dim_loc_ineq;
m           = m_mom + m_loc + m_loc_ineq + 1; 

coef_mom    = coef_all(:,dim_loc_eq+1:dim_loc_eq+n0^2);
coef_mom    = coef_mom';

B           = {};
B_normalize = {};
A           = {};

fprintf('Building B and A... Progress ')
for i = 1:nterms
    if rem(i,10000) == 1
        fprintf('%d/%d ',i,nterms);
    end

    [row,~,~]   = find(coef_mom(:,i));
    SDP_coli    = floor((row-1)./n0) + 1;
    SDP_rowi    = mod(row-1,n0) + 1;
    nnz         = length(SDP_rowi);
    
    Bi          = sparse(SDP_rowi,SDP_coli,ones(nnz,1),n0,n0);
    B{end+1}    = Bi;
    B_normalize{end+1} = Bi/nnz;
    
    mask_triu   = (SDP_rowi >= SDP_coli);
    si          = SDP_rowi(mask_triu);
    sj          = SDP_coli(mask_triu);

    nnz_triu    = length(si);
    
    if nnz_triu > 1
        [~,base_idx]        = max(sj);
        si_base             = si(base_idx);
        sj_base             = sj(base_idx);
        
        si_nonbase          = si; 
        si_nonbase(base_idx)= [];
        sj_nonbase          = sj; 
        sj_nonbase(base_idx)= [];
        
        is_base_diag        = (si_base == sj_base);
        
        if is_base_diag
            A_si            = [si_base];
            A_sj            = [sj_base];
            A_v             = [1];
        else
            A_si            = [si_base,sj_base];
            A_sj            = [sj_base,si_base];
            A_v             = [0.5,0.5];
        end
        
        for nonbase_idx = 1:length(si_nonbase)
            is_nonbase_diag = (si_nonbase(nonbase_idx) == sj_nonbase(nonbase_idx));
            if is_nonbase_diag
                A_sii       = [A_si,si_nonbase(nonbase_idx)];
                A_sjj       = [A_sj,sj_nonbase(nonbase_idx)];
                A_vv        = [A_v,-1];
            else
                A_sii       = [A_si,si_nonbase(nonbase_idx),sj_nonbase(nonbase_idx)];
                A_sjj       = [A_sj,sj_nonbase(nonbase_idx),si_nonbase(nonbase_idx)];
                A_vv        = [A_v,-0.5,-0.5];
            end
            A_temp          = sparse(A_sii,A_sjj,A_vv,n0,n0);
            
            A{end+1}        = A_temp;
        end
    end
end
fprintf('Done.\n')
assert(length(A) == m_mom,'length(A)+length(B) == ndelta!');

%% Now build A's associated with localizing constraints
if dim_loc_eq == 0
    % Do nothing
    A_local = {};
else
    coef_loc    = coef_all(:,1:dim_loc_eq);
    A_local     = {};
    fprintf('Building localizing constraints A_local... Progress ')
    for i = 1:dim_loc_eq
        if rem(i,10000) == 1
            fprintf('%d/%d ',i,m_loc);
        end
        
        [rowi,~,vi] = find(coef_loc(:,i));
        
        Ai      = sparse(n0,n0);
        for j   = 1:length(rowi)
            Ai  = Ai + vi(j) * B_normalize{rowi(j)};
        end
        A_local{end+1} = Ai;
    end
end
fprintf('Done.\n')

%% Leading A
A0 = sparse([1],[1],[1],n0,n0);

%% Combine all A for the main block
A = [{A0},A_local,A];

%% Now build the cost matrix
coef_cost   = coef_all(:,dim_loc_eq+n0^2+1);
[row,~,v]   = find(coef_cost);
C           = sparse(n0,n0);
fprintf('Building cost matrix C... Progress ')
for i = 1:length(row)
    if rem(i,1000) == 1
        fprintf('%d/%d ',i,length(row));
    end
    C       = C + v(i) * B_normalize{row(i)};
end
fprintf('Done.\n')

% Now build the sub PSD constraint g_c (multiple of them)
Ac = {};
for k = 1:length(g_c)
    idx         = dim_loc_eq+n0^2+1+blkIndices(k,n01^2);
    coef_ineq   = coef_all(:,idx); % nterms by n2^2
    A2          = {};
    for ii = 1:n01^2
        row     = mod(ii-1,n01) + 1;
        col     = floor((ii-1)./n01) + 1;
        if row < col
            % Do nothing for upper triangular parts 
        else
            if row == col
                A2tmp       = sparse([row],[col],[-1],n01,n01);
            else
                A2tmp       = sparse([row,col],[col,row],[-0.5,-0.5],n01,n01);
            end
            [termIds,~,v]   = find(coef_ineq(:,ii));
            Atmp            = sparse(n0,n0);
            for iii = 1:length(termIds)
                Atmp        = Atmp + v(iii) * B_normalize{termIds(iii)};
            end
            A{end+1}        = Atmp;
            A2{end+1}       = A2tmp;
        end
    end
    Ac{end+1} = A2;
end
A0 = A;
b0 = sparse(1,1,1,m,1);
C0 = {C};
for k = 1:length(g_c)
    C0 = [C0;{sparse(n01,n01)}];
end


%% Now build the chordal blocks
%% the 1-N blocks v = [1;x;c;cr;theta(i);theta(i)*x;theta(i)*cr]
%% Since there is an inequality constraint too, it will generate 2*N blocks
Aall            = {};
A1all           = {};
ball            = [];
Call            = {};
A0append        = {};
for blkidx = 1:N
    basis       = [1;x;c;cr;theta(blkidx);theta(blkidx)*x;theta(blkidx)*cr];
    n           = length(basis);
    basis_r     = get_multiplier_basis([x;c;theta(blkidx)],basis,h_r(1));
    basis_theta = get_multiplier_basis([x;c;theta(blkidx)],basis,h_theta(blkidx));
    basis_g     = [1;theta(blkidx)];
    
    out         = gen_chordal_subblk_catreg(...
                    basis,basis_r,h_r,basis_theta,h_theta(blkidx),g_t,basis_g,f_cost(blkidx+1));
    Acell       = out.A;
    A1          = out.A1;
    b           = out.b;
    C           = out.C;
    
    % add constraint that the top-left block is the same as
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
    A1all{end+1}     = A1;
end

%% Convert to standard SDPT3 format
b           = [b0;ball];
blk         = cell(2*N+1+length(g_c),2);
blk{1,1}    = 's'; blk{1,2} = n0;
for k = 1:length(g_c)
    blk{k+1,1} = 's'; blk{k+1,2} = n01;
end
n1          = out.blk{2,2};
n1delta     = triangle_number(n1);
ndelta      = triangle_number(n);
for i = 1:N 
    blk{2*i+length(g_c),1} = 's';
    blk{2*i+length(g_c),2} = n;
    blk{2*i+length(g_c)+1,1} = 's';
    blk{2*i+length(g_c)+1,2} = n1;
end

A0t     = sparsesvec(blk(1,:),A0);
for i = 1:N
    A0t = [A0t,...
           sparse(n0delta,out.m),...
           sparsesvec(blk(1,:),A0append{i})];
end
Act     = {};
for k = 1:length(g_c)
    Actk = [sparse(n01delta,m - m_loc_ineq),...
            sparse(n01delta,(k-1)*n01delta),...
            sparsesvec(blk(k+1,:),Ac{k}),...
            sparse(n01delta,(length(g_c)-k)*n01delta),...
            sparse(n01delta,N*(out.m + n0delta))];
    Act  = [Act;{Actk}];
end
At      = [{A0t};Act];

for i = 1:N
    Ait = [sparse(ndelta,length(b0)),...
           sparse(ndelta,(i-1)*length(Aall{i})),...
           sparsesvec(blk(2*i+length(g_c),:),Aall{i}),...
           sparse(ndelta,(N-i)*length(Aall{i}))];
       
    A1it = [sparse(n1delta,length(b0)),...
            sparse(n1delta,(i-1)*length(Aall{i})),...
            sparse(n1delta,out.m_mom+out.m_loc),...
            sparsesvec(blk(2*i+length(g_c)+1,:),A1all{i}),...
            sparse(n1delta,n0delta),...
            sparse(n1delta,(N-i)*length(Aall{i}))];
    At  = [At;{Ait};{A1it}];
end

SDP.blk = blk;
SDP.At  = At;
SDP.m   = length(b);
SDP.C   = [C0;Call];
SDP.b   = b;

tf    = toc(t0);
fprintf('\nDone in %g seconds.\n',tf);
fprintf('===================================================================\n')

end