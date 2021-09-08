function [Acell,b,C] = poly2SDP(n,degmat,coef_all,dim_loc,varargin)

params = inputParser;
params.CaseSensitive = false;
params.addParameter('verbose',1, @(x) isscalar(x));
params.parse(varargin{:});

verbose     = params.Results.verbose;

ndelta      = triangle_number(n);
% there are nterms unique monomials, which means we will have
% nterms B matrices, where each Bi indicate the location of 
% the i-th monomial in the moment matrix v*v';
% This also means the moment matrix p*p' is a linear subspace of
% dimension nterms in sym{n}
nterms      = size(degmat,1);
% Because p*p' is a linear subspace of dimension nterms,
% it's orthogonal complement will be a linear subspace of
% dimension ndelta - nterms
m_mom       = ndelta - nterms;
% Total number of linear constraints
m_loc       = dim_loc;
m           = m_mom + m_loc + 1; 

if verbose > 0
    fprintf('SDP: n = %d, m = %d, m_mom = %d, m_loc = %d, ndelta = %d.\n',n,m,m_mom,m_loc,ndelta);
end

coef_mom    = coef_all(:,dim_loc+1:end-1);
coef_mom    = coef_mom';

B           = {};
B_normalize = {};
A           = {};

if verbose > 0
    fprintf('Building B and A... Progress ')
end
for i = 1:nterms
    if verbose > 0
        if rem(i,10000) == 1
            fprintf('%d/%d ',i,nterms);
        end
    end

    [row,~,~]   = find(coef_mom(:,i));
    SDP_coli    = floor((row-1)./n) + 1;
    SDP_rowi     = mod(row-1,n) + 1;
    nnz         = length(SDP_rowi);
    
    Bi          = sparse(SDP_rowi,SDP_coli,ones(nnz,1),n,n);
    % assert(norm(Bi-Bi','fro')==0,'B should be symmetric!');
    B{end+1}    = Bi;
    B_normalize{end+1} = Bi/nnz;
    
    mask_triu   = (SDP_rowi >= SDP_coli);
    si          = SDP_rowi(mask_triu);
    sj          = SDP_coli(mask_triu);

    nnz_triu    = length(si);
    
    if nnz_triu > 1
        [~,base_idx]     = max(sj);
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
            % assert(norm(A_temp-A_temp','fro')==0,'A should be symmetric!');
            A{end+1}        = A_temp;
        end
    end
end
if verbose > 0
    fprintf('Done.\n')
end
assert(length(A)==m_mom,'length(A)+length(B) == ndelta!');

%% Now build A's associated with localizing constraints
if dim_loc == 0
    % Do nothing
    A_local = {};
else
coef_loc    = coef_all(:,1:dim_loc);
A_local     = {};

if verbose > 0
    fprintf('Building localizing constraints A_local... Progress ')
end
for i = 1:dim_loc
    if verbose > 0
        if rem(i,10000) == 1
            fprintf('%d/%d ',i,m_loc);
        end
    end
    
    [rowi,~,vi] = find(coef_loc(:,i));
    
    Ai      = sparse(n,n);
    for j   = 1:length(rowi)
        Ai  = Ai + vi(j) * B_normalize{rowi(j)};
    end
    A_local{end+1} = Ai;
end
end
if verbose > 0
    fprintf('Done.\n')
end

%% Leading A
A0 = sparse([1],[1],[1],n,n);

%% Combine all A
Acell = [{A0},A_local,A];
b = sparse([1],[1],[1],length(Acell),1);
assert(m == length(Acell), 'SDP equality number wrong.');

%% Now build the cost matrix
coef_cost   = coef_all(:,end);
[row,~,v] = find(coef_cost);
C           = sparse(n,n);
% fprintf('Building cost matrix C... Progress ')
for i = 1:length(row)
    % if rem(i,1000) == 1
    %     fprintf('%d/%d ',i,length(row));
    % end
    C       = C + v(i) * B_normalize{row(i)};
end
% fprintf('Done.\n')

end