function SDP = relax_absolute_pose_estimation_v1(problem,varargin)
%% Apply a sparse second-order relaxation to absolute pose estimation
%% Full perspective camera, R and t relaxed together
%% Depending on multivariate polynomial package in SPOT
%% Translation bounded modeled as an extra PSD constraint
%% t(3) >= 0 modelled as an extra PSD constraint
%% tan(FOV/2)^2 * t(3)^2 >= (t(1)^2 + t(2)^2); in FOV CONE (July 01)
%% Heng Yang, July 01, 2021

params = inputParser;
params.CaseSensitive = false;
params.addParameter('checkMonomials',true, @(x) islogical(x));
params.parse(varargin{:});

checkMonomials = params.Results.checkMonomials;

fprintf('\n===================================================================')
fprintf('\nApplying SDP relaxation to absolute pose estimation')
fprintf('\n===================================================================\n')
t0              = tic;

N               = problem.N;
meas3D          = problem.X;
meas2D          = problem.x;
noiseBoundSq    = problem.noiseBoundSq;
tBound          = problem.translationBound;
tBoundSq        = tBound^2; % t'*t <= tBoundSq
dBound          = problem.depthBound;
barc2           = 1.0;
FOV             = deg2rad(problem.FOV);

%% define POP variables
nrPrimalVars    = 9+3 + N; % rotation: 9, translation: 3, binary: N
p               = msspoly('p',nrPrimalVars);
r               = p(1:9);
R               = reshape(r,3,3); 
col1            = r(1:3);
col2            = r(4:6);
col3            = r(7:9);
row1            = R(1,:)';
row2            = R(2,:)';
row3            = R(3,:)';
t               = p(10:12);
x               = [r;t];
theta           = p(13:nrPrimalVars);

%% compute the cost function
residuals = {};
for i = 1:N 
    bearingi    = [meas2D(:,i);1];
    bearingi    = bearingi / norm(bearingi);
    pointi      = R * meas3D(:,i) + t;
    res         = pointi' * (eye(3) - bearingi*bearingi') * pointi; % point to line distance
    residuals{end+1} ...
                = res / noiseBoundSq;
end
f_cost = 0;
for i = 1:N 
    f_cost = f_cost + (1+theta(i))/2 * residuals{i} + (1-theta(i))/2 * barc2;
end

%% define equality constraints
h_x  = [1.0-col1'*col1;...
       1.0-col2'*col2;...
       1.0-col3'*col3;... % columns unit length
       1.0-row1'*row1;...
       1.0-row2'*row2;...
       1.0-row3'*row3;... % rows unit length
       col1'*col2;...
       col2'*col3;...
       col3'*col1;... % columns orthogonal
       row1'*row2;...
       row2'*row3;...
       row3'*row1;... % rows orthogonal
       cross(row1,row2) - row3;...
       cross(row2,row3) - row1;...
       cross(row3,row1) - row2;... % rows righthandedness
       cross(col1,col2) - col3;...
       cross(col2,col3) - col1;...
       cross(col3,col1) - col2]; % columns righthandedness

% h_x  = [1.0-col1'*col1;...
%         1.0-col2'*col2;...
%         1.0-col3'*col3;... % column unit length
%         col1'*col2;...
%         col2'*col3;...
%         col3'*col1;... % colums orthogonal
%         cross(col1,col2) - col3;...
%         cross(col2,col3) - col1;...
%         cross(col3,col1) - col2]; % columns righthandedness

h_theta = [];
for i = 1:N 
    h_theta =[h_theta; 1-theta(i)^2];
end

g_x = [tBoundSq - t'*t; ... % Translation bounded
       t(3)-dBound; ...  % in front of camera
       (tan(FOV/2))^2 * t(3)^2 - (t(1)^2 + t(2)^2)]; % in FOV cone

%% Formulate the sparse second-order relaxation
basis_p     = [1;x;theta;mykron(theta,x)];
n           = length(basis_p);
basis_x     = monomials(theta,0:2);
basis_theta = monomials(x,0:2);
basis_g     = [1;theta];
n1          = length(basis_g);
fprintf('Computing localizing and moment polynomials ...')
time_start  = tic;
tmp = mykron(basis_g,basis_g);
pop = [mykron(basis_x,h_x);...
       mykron(basis_theta,h_theta);...
       mykron(basis_p,basis_p);...
       f_cost;...
       g_x(1)*tmp;...
       g_x(2)*tmp;...
       g_x(3)*tmp];
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
n1delta     = triangle_number(n1);
m_loc_ineq  = length(g_x)*n1delta; % each inequality constraint leads to n1delta equality constraints
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
coef_ineq   = coef_all(:,dim_loc_eq+n^2+1+1:dim_loc_eq+n^2+1+n1^2); % nterms by n1^2
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

%% Now build the constraints from the third PSD constraint
coef_ineq   = coef_all(:,dim_loc_eq+n^2+1+n1^2+1:dim_loc_eq+n^2+1+2*n1^2); % nterms by n1^2
A2          = {};
for ii = 1:n1^2
    row     = mod(ii-1,n1) + 1;
    col     = floor((ii-1)./n1) + 1;
    if row < col
        % Do nothing for upper triangular parts 
    else
        if row == col
            A2tmp       = sparse([row],[col],[-1],n1,n1);
        else
            A2tmp       = sparse([row,col],[col,row],[-0.5,-0.5],n1,n1);
        end
        [termIds,~,v]   = find(coef_ineq(:,ii));
        Atmp            = sparse(n,n);
        for iii = 1:length(termIds)
            Atmp        = Atmp + v(iii) * B_normalize{termIds(iii)};
        end
        A{end+1}        = Atmp;
        A2{end+1}       = A2tmp;
    end
end

%% Now build the constraints from the fourth PSD constraint
coef_ineq   = coef_all(:,dim_loc_eq+n^2+1+2*n1^2+1:dim_loc_eq+n^2+1+3*n1^2); % nterms by n1^2
A3          = {};
for ii = 1:n1^2
    row     = mod(ii-1,n1) + 1;
    col     = floor((ii-1)./n1) + 1;
    if row < col
        % Do nothing for upper triangular parts 
    else
        if row == col
            A3tmp       = sparse([row],[col],[-1],n1,n1);
        else
            A3tmp       = sparse([row,col],[col,row],[-0.5,-0.5],n1,n1);
        end
        [termIds,~,v]   = find(coef_ineq(:,ii));
        Atmp            = sparse(n,n);
        for iii = 1:length(termIds)
            Atmp        = Atmp + v(iii) * B_normalize{termIds(iii)};
        end
        A{end+1}        = Atmp;
        A3{end+1}       = A3tmp;
    end
end

assert(length(A) == m,'Total number of equality constraints wrong.')
assert(length(A1)+length(A2)+length(A3)== m_loc_ineq,'Equality constraints from sub PSD blks wrong.')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
blk{1,1}        = 's';
blk{1,2}        = n;
blk{2,1}        = 's';
blk{2,2}        = n1;
blk{3,1}        = 's';
blk{3,2}        = n1;
blk{4,1}        = 's';
blk{4,2}        = n1;
b               = sparse([1],[1],[1],m,1);
%% svec in sdpt3 format and output standard data
At0             = sparsesvec(blk(1,:),A);

At1             = [sparse(n1delta,m-m_loc_ineq),...
                   sparsesvec(blk(2,:),A1),...
                   sparse(n1delta,n1delta),...
                   sparse(n1delta,n1delta)];
               
At2             = [sparse(n1delta,m-m_loc_ineq),...
                   sparse(n1delta,n1delta),...
                   sparsesvec(blk(3,:),A2),...
                   sparse(n1delta,n1delta)];
               
At3             = [sparse(n1delta,m-m_loc_ineq),...
                   sparse(n1delta,n1delta),...
                   sparse(n1delta,n1delta),...
                   sparsesvec(blk(4,:),A3),];
               
At              = {At0;At1;At2;At3};

SDP.relaxversion = 'v1';
SDP.blk = blk;
SDP.At  = At;
SDP.m   = m;
SDP.C   = {C;sparse(n1,n1);sparse(n1,n1);sparse(n1,n1)};
SDP.b   = b;
% SDP.M   = (4+tBoundSq)*(N+1); % Bound on the trace of the moment matrix

SDP.M   = [(4+tBoundSq)*(N+1);...
           tBoundSq * (1+N);...
           sqrt(tBoundSq) * (1+N);...
           (tan(FOV/2))^2 * tBoundSq * (1+N)];

tf    = toc(t0);
fprintf('\nDone in %g seconds.\n',tf);
fprintf('===================================================================\n')
end