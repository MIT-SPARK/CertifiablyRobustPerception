function [SDP,info] = dense_sdp_relax(problem,kappa_user)
%% Implements dense Lasserre's hierarchy
%% Heng Yang
%% Sept. 19, 2021
%% Oct. 13, 2021: fixed a bug in generating redundant equality constraints

if nargin < 2; kappa_user = 0; end
if ~isfield(problem,'vars'); error('Please provide variables of the POP.'); end
if ~isfield(problem,'objective'); error('Please provide objective function.'); end
if ~isfield(problem,'equality'); problem.equality = []; end
if ~isfield(problem,'inequality'); problem.inequality = []; end

fprintf('\n======================= Dense SDP Relaxation =======================\n')

%% copy POP data and decide relaxation order
x       = problem.vars;
f       = problem.objective;
h       = problem.equality;
g       = problem.inequality;

degf    = deg(f);
if isempty(h); degh = 0; else; degh = deg(h); end
if isempty(g); degg = 0; else; degg = deg(g); end
maxdeg  = max([degf;degh;degg]);

kappa_min   = ceil(maxdeg/2);
kappa       = max(kappa_min,kappa_user);

if kappa == 0; error('kappa=0, looks like this is not a valid POP.'); end
if kappa == 1; v = [1;x]; end
if kappa > 1; v = [1;x;monomials(x,2:kappa)]; end

J = diff(f,x);
info.f                   = f;
info.J                   = J'; % for evaluate cost gradient
info.v                   = v; % for lifting
info.d                   = length(x);

%% compute multipliers
fprintf('Compute equality and inequality multipliers ... ')
kappa2  = 2*kappa;
pop     = [mykron(v,v);f];
% equalities
l_h     = length(h);
dim_h   = 0;
if l_h > 0
    for i = 1:l_h
        hi      = h(i);
        deghi   = deg(hi);
        
        if deghi == 0; error('One of your equality constraints has zero degree.'); end
        
        lamhi   = monomials(x,0:(kappa2 - deghi));
        pop     = [pop;hi*lamhi];
        dim_h   = dim_h + length(lamhi);
    end
end
% inequalities
l_g     = length(g);
lam_g   = {};
n1s     = [];
if l_g > 0 % if there are inequality constraints
    for i = 1:l_g
        gi      = g(i);
        deggi   = deg(gi);
        
        if deggi == 0; error('One of your inequality constraints has zero degree.'); end
        
        ordergi = kappa - ceil(deggi/2); % order of the multiplier moment matrix
        
        if ordergi == 0
            lamgi = 1;
        elseif ordergi == 1
            lamgi = [1;x];
        else
            lamgi = [1;x;monomials(x,2:ordergi)];
        end
        
        lamgi_flat  = mykron(lamgi,lamgi);
        pop         = [pop;gi*lamgi_flat];
        lam_g       = [lam_g;{lamgi}];
        n1s         = [n1s;length(lamgi)];
    end
end
fprintf('Done.\n')

fprintf('POP | maxdeg: %d, l_h: %d, l_g: %d, kappa: %d.\n',...
    maxdeg,l_h,l_g,kappa);

%% decompose polynomials to generate SDP data
[~,degmat,coef_all] = decomp(pop);
coef_all            = coef_all';

n                   = length(v); % size of the moment matrix
ndelta              = triangle_number(n);
nterms              = size(degmat,1);
m_mom               = ndelta - nterms;
m_h                 = dim_h; % number of constraints due to equality
m_g                 = sum( triangle_number(n1s) );
m                   = m_mom + m_h + m_g + 1;

fprintf('SDP | blknum: %d, moment mat: %d, m: %d, m_mom: %d, m_h: %d, m_g: %d.\n',...
    length(n1s)+1,n,m_mom+m_h+m_g+1,m_mom,m_h,m_g);

coef_mom    = coef_all(:,1:n^2);
coef_mom    = coef_mom';

B           = {};
B_normalize = {};
A           = {};

fprintf('Build moment constraints ... ')
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

%% Build A's associated with localizing equality constraints
if m_h == 0
    % Do nothing
    A_local = {};
else
    coef_loc    = coef_all(:,1+n^2+(1:m_h));
    A_local     = {};
    fprintf('Build localizing equality constraints ... ')
    for i = 1:m_h
        if rem(i,10000) == 1
            fprintf('%d/%d ',i,m_h);
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
% Combine all A for the main block
A = [{A0},A_local,A];

%% Now build the cost matrix
coef_cost   = coef_all(:,n^2+1);
[row,~,v]   = find(coef_cost);
C           = sparse(n,n);
fprintf('Build cost matrix C ... ')
for i = 1:length(row)
    if rem(i,1000) == 1
        fprintf('%d/%d ',i,length(row));
    end
    C       = C + v(i) * B_normalize{row(i)};
end
fprintf('Done.\n')

%% Now build the localizing inequality constraints
Ag  = {};
shift = 0;
if l_g > 0
    fprintf('Build localizing inequality constraints ... ')
    for k = 1:l_g
        n1          = n1s(k);
        coef_ineq   = coef_all(:,1+n^2+m_h+shift+(1:n1^2)); % nterms by n1^2
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
        Ag{end+1} = A1;
        shift     = shift + n1^2;
    end
    fprintf('Done.\n')
end

assert(length(A) == m,'Total number of equality constraints wrong.')

%% convert to SDPT3 format
fprintf('Generate SDPT3 data ... ')
blk{1,1}        = 's';
blk{1,2}        = n;
for k = 1:l_g
    blk{1+k,1}  = 's';
    blk{1+k,2}  = n1s(k);
end
b               = sparse([1],[1],[1],m,1);
At0             = sparsesvec(blk(1,:),A);
At              = {At0};
n1deltas        = triangle_number(n1s);
for k = 1:l_g
    n1delta     = triangle_number(n1s(k));
    Atg         = [sparse(n1delta,m-m_g),... % moment and localizing equality
                   sparse(n1delta,sum(n1deltas(1:k-1))),... % sub PSD g
                   sparsesvec(blk(1+k,:),Ag{k}),... % sub PSD k-th g
                   sparse(n1delta,sum(n1deltas(k+1:end)))]; 
    At          = [At;{Atg}];
end
C       = {C};
for k = 1:l_g
    C   = [C;{sparse(n1s(k),n1s(k))}];
end

SDP.blk = blk;
SDP.At  = At;
SDP.C   = C;
SDP.b   = b;

fprintf('Done.\n')

%% convert to sedumi format
fprintf('Generate sedumi data ... ')
sK.s  = [n;n1s];

A0t     = sparsevec(blk(1,:),A);
At      = {A0t};
for k = 1:l_g
    n1sq        = (n1s(k))^2;
    Atg         = [sparse(n1sq,m-m_g),... % moment and localizing
                   sparse(n1sq,sum(n1deltas(1:k-1))),... % sub PSD g
                   sparsevec(blk(1+k,:),Ag{k}),... % sub PSD k-th g
                   sparse(n1sq,sum(n1deltas(k+1:end)))]; 
    At          = [At;{Atg}];
end

sdata.K     = sK;
sdata.At    = cat(1,At{:});
sdata.b     = b;

sc          = [];
for i = 1:length(SDP.C)
    sc      = [sc;sparsevec(blk(i,:),SDP.C(i))];
end
sdata.c     = sc;

SDP.sedumi   = sdata;

fprintf('Done.\n')
fprintf('====================================================================\n\n\n')
info.kappa  = kappa;
info.lam_g  = lam_g;

end