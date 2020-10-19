function projection = get_dual_projection_single_rotation_averaging(N)
    % get projection maps
    t0 = tic;
    nrPrimalVars = 9 + N;
    p = mpvar('p',[nrPrimalVars,1]);
    r = p(1:9); R = reshape(r,[3,3]);
    r1 = R(:,1); r2 = R(:,2); r3 = R(:,3);
    theta = p(10:end);
    % build my constraints
    g_r = [1-r1'*r1;1-r2'*r2;1-r3'*r3;...
           r1'*r2;r2'*r3;r3'*r1;...
           hatmap(r1)*r2-r3;hatmap(r2)*r3-r1;hatmap(r3)*r1-r2];
    g_theta = [];
    for i = 1:N 
        g_theta =[g_theta; 1-theta(i)^2];
    end
    nrCons_r = length(g_r);
    nrCons_theta = length(g_theta);

    % use d to denote dual variables (PSD matrices, lambdas)
    % f(p) - f_est - \sum lambda_ri*g_ri - \sum lambda_thetai*g_thetai = (basis_p' * Xp * basis_p)
    % lambda_ri = d_ri' * basis_r
    % lambda_thetai = d_thetai' * basis_theta
    basis_r = monomials(theta,0:2);
    basis_theta = monomials(r,0:2);
    basis_p = [ monomials(p,0:1); kron(theta,r) ];
    dim_Xp = length(basis_p);
    
    nrDualVars_r = length(basis_r) * nrCons_r;
    nrDualVars_theta = length(basis_theta) * nrCons_theta;
    nrDualVars_p = dim_Xp^2;
    nrDualVars_free = nrDualVars_r + nrDualVars_theta;
    nrDualVars_psd = nrDualVars_p;
    nrDualVars = nrDualVars_free + nrDualVars_psd;
    % Build the affine subspace [A1, A2] * d = b
    poly_fuse = [kron(g_r,basis_r);...
                 kron(g_theta,basis_theta);...
                 kron(basis_p,basis_p)];

    assert(length(poly_fuse) == nrDualVars, 'Dimension is wrong.');

    A = poly_fuse.coefficient;
    A1 = A(:,1:nrDualVars_free);
    A2_raw = A(:,nrDualVars_free+1:end)';
    nrTerms = size(A,1);
    nrTerms_full = nchoosek(nrPrimalVars + 2*2, 2*2);
    % convert to svec format to speed up
    dim_svecp = dim_Xp * (dim_Xp + 1)/2;
    nrDualVars_svec = nrDualVars_free + dim_svecp;
    A2 = sparse(dim_svecp,nrTerms);
    for i = 1:size(A2_raw,2)
        A2(:,i) = svec(reshape(A2_raw(:,i),[dim_Xp,dim_Xp]));
    end
    A2 = A2';
    A = [A1, A2];

    [inverse_A, inverse_preb] = get_dual_projection(A1,A2);
    projection.type = 'single rotation averaging';
    projection.N = N;
    projection.A = A;
    projection.inv_A = inverse_A;
    projection.pre_b = inverse_preb;
    projection.degmat = poly_fuse.degmat;
    projection.varname = poly_fuse.varname;
    projection.nrDualVars_free = nrDualVars_free;
    projection.dim_svecp = dim_svecp;
    projection.dim_Xp = dim_Xp;

    t_total = toc(t0);

    fprintf('projection single rotation averaging: N=%d.\n',N);
    fprintf('# of Terms w/ reduction: \t %d\n',nrTerms);
    fprintf('# of Terms w/o reduction: \t %d\n',nrTerms_full);
    fprintf('# of Columns in A1: \t %d\n',size(A1,2));
    fprintf('# of Columns in A2: \t %d (Orthogonal)\n',size(A2,2));
    fprintf('Obtained projection map in %g seconds.\n',t_total);

    % save the projection map
    % assume this function is called from certifiable-perception/example_**
    % directories
    filename = sprintf('../solver/single_rotation_averaging/projection_maps/projection_sra_%d.mat',N);
    save(filename,'projection');
end

function [proj_dual_A,proj_dual_pre_b] = get_dual_projection(A1,A2)
    [nrTerms,~] = size(A1);
    A = [A1,A2];
    P = A2*A2';
    P_inv_diag = 1./diag(P);
    P_inv = spdiags(P_inv_diag,0,nrTerms,nrTerms);
    Z = speye(size(A1,2)) + A1'*P_inv*A1;
    M_inv = P_inv - P_inv * A1 * (Z \ A1') * P_inv;

    proj_dual_pre_b = A'*M_inv;
    proj_dual_A = speye(size(A,2)) - proj_dual_pre_b*A;
end