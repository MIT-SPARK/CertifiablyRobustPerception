function SDP = relax_mesh_registration(problem,varargin)
%% Apply a sparse second-order relaxation to mesh registration
%% Depending on multivariate polynomial package
%% Translation bounded modeled as an extra PSD constraint

params = inputParser;
params.CaseSensitive = false;

params.addParameter('checkMonomials',true,...
    @(x) islogical(x));
params.addParameter('tighterRelax',false,...
    @(x) islogical(x));

params.parse(varargin{:});

checkMonomials = params.Results.checkMonomials;
tighterRelax   = params.Results.tighterRelax;

fprintf('\n===================================================================')
fprintf('\nApplying SDP relaxation to mesh registration problem')
fprintf('\n===================================================================\n')
t0              = tic;

N               = problem.N;
normalM         = problem.normalM;
normalP         = problem.normalP;
pointM          = problem.pointM;
pointP          = problem.pointP;
pointNoiseBoundSq  = problem.pointNoiseBoundSq;
normalNoiseBoundSq = problem.normalNoiseBoundSq;
tBound          = problem.translationBound;
tBoundSq        = tBound^2; % t'*t <= tBoundSq
barc2           = 1.0;

%% define POP variables
nrPrimalVars    = 9+3 + N; % rotation: 9, translation: 3, binary: N
p               = msspoly('p',nrPrimalVars);
r               = p(1:9);
R               = reshape(r,3,3); 
col1            = r(1:3);
col2            = r(4:6);
col3            = r(7:9);
t               = p(10:12);
x               = [r;t];
theta           = p(13:nrPrimalVars);

%% compute the cost function
residuals = {};
for i = 1:N
    pointM_i    = pointM(:,i);
    pointP_i    = pointP(:,i);
    normalM_i   = normalM(:,i);
    normalP_i   = normalP(:,i);
    residual_point_i  = ( normalM_i' * (R*pointP_i - pointM_i + t) )^2;
    residual_normal_i = normalP_i'*normalP_i + normalM_i'*normalM_i - 2*normalP_i'*R'*normalM_i;
    residuals{end+1} ...
                = (residual_point_i/pointNoiseBoundSq + residual_normal_i/normalNoiseBoundSq)/2;
end
f_cost = 0;
for i = 1:N 
    f_cost = f_cost + (1+theta(i))/2 * residuals{i} + (1-theta(i))/2 * barc2;
end

%% define equality constraints
h_x = [1.0-col1'*col1;...
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

g_x = tBoundSq - t'*t; % Translation bounded

%% Formulate the sparse second-order relaxation
if tighterRelax
    basis_p = [1;x;theta;monomials(x,2);mykron(theta,x)];
else % default
    basis_p = [1;x;theta;mykron(theta,x)];
end
n       = length(basis_p);
if tighterRelax
    basis_x     = monomials(p,0:2);
else
    basis_x     = monomials(theta,0:2);
end
basis_theta = monomials(x,0:2);
if tighterRelax
    basis_g     = [1;p];
else
    basis_g     = [1;theta];
end
n1          = length(basis_g);

fprintf('Computing localizing and moment polynomials ...')
time_start  = tic;
pop = [mykron(basis_x,h_x);...
       mykron(basis_theta,h_theta);...
       mykron(basis_p,basis_p);
       f_cost;
       g_x*mykron(basis_g,basis_g)];
[~,degmat,coef_all] = decomp(pop);
coef_all            = coef_all';
time_prep   = toc(time_start);
fprintf(' Done in %g seconds.\n',time_prep);

if checkMonomials
    fprintf('Checking consistency of monomials ...')
    time_check0   = tic;
    monomials_mom = mono(mykron(basis_p,basis_p));
    fprintf(' %d ... %d ...',size(degmat,1),length(monomials_mom));
    assert(size(degmat,1) == length(monomials_mom),'monomials not consistent');
    time_check    = toc(time_check0);
    fprintf('Done in %g seconds.\n',time_check);
end

dim_loc_eq = length(basis_x)*length(h_x) + length(basis_theta)*length(h_theta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generate standard SDP data from degmat and coefficients
ndelta      = triangle_number(n);
nterms      = size(degmat,1);
m_mom       = ndelta - nterms;
m_loc       = dim_loc_eq;
m_loc_ineq  = triangle_number(n1); % The second PSD constraint leads to triangle_number(n1) linear equality constraints
m           = m_mom + m_loc + m_loc_ineq + 1; 

fprintf('SDP: n = %d, n1 = %d, m = %d, m_mom = %d, m_loc = %d, m_loc_ineq = %d, ndelta = %d.\n',...
        n,n1,m,m_mom,m_loc,m_loc_ineq,ndelta);

coef_mom    = coef_all(:,dim_loc_eq+1:dim_loc_eq+n^2);
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
    SDP_coli    = floor((row-1)./n) + 1;
    SDP_rowi    = mod(row-1,n) + 1;
    nnz         = length(SDP_rowi);
    
    Bi          = sparse(SDP_rowi,SDP_coli,ones(nnz,1),n,n);
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
            A_temp          = sparse(A_sii,A_sjj,A_vv,n,n);
            
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
        
        Ai      = sparse(n,n);
        for j   = 1:length(rowi)
            Ai  = Ai + vi(j) * B_normalize{rowi(j)};
        end
        A_local{end+1} = Ai;
    end
end
fprintf('Done.\n')

%% Leading A
A0 = sparse([1],[1],[1],n,n);

%% Combine all A for the main block
A = [{A0},A_local,A];

%% Now build the cost matrix
coef_cost   = coef_all(:,dim_loc_eq+n^2+1);
[row,~,v]   = find(coef_cost);
C           = sparse(n,n);
fprintf('Building cost matrix C... Progress ')
for i = 1:length(row)
    if rem(i,1000) == 1
        fprintf('%d/%d ',i,length(row));
    end
    C       = C + v(i) * B_normalize{row(i)};
end
fprintf('Done.\n')

%% Now build the constraints from the second PSD constraint
coef_ineq   = coef_all(:,dim_loc_eq+n^2+1+1:end); % nterms by n1^2
A1          = {};
for ii = 1:n1^2
    row     = mod(ii-1,n1) + 1;
    col     = floor((ii-1)./n1) + 1;
    if row < col
        % Do nothing for upper triangular parts 
    else
        if row == col
            A1tmp       = sparse([row],[col],[-1],n1,n1);
        else
            A1tmp       = sparse([row,col],[col,row],[-0.5,-0.5],n1,n1);
        end
        [termIds,~,v]   = find(coef_ineq(:,ii));
        Atmp            = sparse(n,n);
        for iii = 1:length(termIds)
            Atmp        = Atmp + v(iii) * B_normalize{termIds(iii)};
        end
        A{end+1}        = Atmp;
        A1{end+1}       = A1tmp;
    end
end

assert(length(A) == m,'Total number of equality constraints wrong.')
assert(length(A1)== m_loc_ineq,'Equality constraints from second PSD blk wrong')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blk{1,1}        = 's';
blk{1,2}        = n;
blk{2,1}        = 's';
blk{2,2}        = n1;
b               = sparse([1],[1],[1],m,1);
%% svec in sdpt3 format and output standard data
At0             = sparsesvec(blk(1,:),A);
At1             = [sparse(triangle_number(n1),m-m_loc_ineq),sparsesvec(blk(2,:),A1)];
At              = {At0;At1};

SDP.blk = blk;
SDP.At  = At;
SDP.m   = m;
SDP.C   = {C;sparse(n1,n1)};
SDP.b   = b;
% SDP.M   = (4+tBoundSq)*(N+1); % bound on the trace of the moment matrix

SDP.M   = [(4+tBoundSq)*(N+1); tBoundSq * (1+N)];

tf    = toc(t0);
fprintf('\nDone in %g seconds.\n',tf);
fprintf('===================================================================\n')
end