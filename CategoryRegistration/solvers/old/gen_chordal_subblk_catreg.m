function out = gen_chordal_subblk_catreg(basis,basis_x,h_x,basis_theta,h_theta,g_x,basis_g,f_cost)

    n       = length(basis);
    n1      = length(basis_g);
    
    
    basis_g_flat = mykron(basis_g,basis_g);
    
    pop          = [mykron(basis_x,h_x);...
                    mykron(basis_theta,h_theta);...
                    mykron(basis,basis);...
                    f_cost;...
                    g_x*basis_g_flat];
    
                
    [pp,degmat,coef_all] = decomp(pop);
    coef_all            = coef_all';
    dim_loc_eq          = length(basis_x)*length(h_x) + length(basis_theta)*length(h_theta);
    
    %% check monomial consistency
    monomials_mom = mono(mykron(basis,basis));
    fprintf(' %d ... %d ...',size(degmat,1),length(monomials_mom));
    assert(size(degmat,1) == length(monomials_mom),'monomials not consistent');
    
    %% Build SDP data
    ndelta      = triangle_number(n);
    nterms      = size(degmat,1);
    m_mom       = ndelta - nterms;
    m_loc       = dim_loc_eq;
    n1delta     = triangle_number(n1);
    m_loc_ineq  = n1delta; % The second PSD constraint leads to triangle_number(n1) linear equality constraints
    m           = m_mom + m_loc + m_loc_ineq; 
        
    coef_mom    = coef_all(:,dim_loc_eq+1:dim_loc_eq+n^2);
    coef_mom    = coef_mom';
    
    B           = {};
    B_normalize = {};
    A           = {};
    
    for i = 1:nterms
    
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
                if degmat(i,end) > 0
                    A{end+1}        = A_temp;
                end
            end
        end
    end
    % assert(length(A) == m_mom,'length(A)+length(B) == ndelta!');
    m_mom   = length(A);
    m       = m_mom + m_loc + m_loc_ineq; 
    fprintf('SDP: n = %d, n1 = %d, m = %d, m_mom = %d, m_loc = %d, m_loc_ineq = %d, ndelta = %d.\n',...
            n,n1,m,m_mom,m_loc,m_loc_ineq,ndelta);
    
    % Now build A's associated with localizing constraints
    if dim_loc_eq == 0
        % Do nothing
        A_local = {};
    else
        coef_loc    = coef_all(:,1:dim_loc_eq);
        A_local     = {};
        
        for i = 1:dim_loc_eq
            [rowi,~,vi] = find(coef_loc(:,i));
            
            Ai      = sparse(n,n);
            for j   = 1:length(rowi)
                Ai  = Ai + vi(j) * B_normalize{rowi(j)};
            end
            A_local{end+1} = Ai;
        end
    end
    
    % Combine all A for the main block
    A = [A_local,A]; A = A';
    
    % Now build the cost matrix
    coef_cost   = coef_all(:,dim_loc_eq+n^2+1);
    [row,~,v]   = find(coef_cost);
    C           = sparse(n,n);
    for i = 1:length(row)
        C       = C + v(i) * B_normalize{row(i)};
    end
    
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
    
    blk{1,1}        = 's';
    blk{1,2}        = n;
    blk{2,1}        = 's';
    blk{2,2}        = n1;
    b               = sparse(m,1);
    
    out.A           = A;
    out.A1          = A1';
    out.blk         = blk;
    out.b           = b;
    out.C           = {C;sparse(n1,n1)};
    out.m           = m;
    out.m_mom       = m_mom;
    out.m_loc       = m_loc;
    out.m_loc_ineq  = m_loc_ineq;
    
    end