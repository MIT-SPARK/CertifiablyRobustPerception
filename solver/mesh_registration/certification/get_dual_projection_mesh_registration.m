function projection = get_dual_projection_mesh_registration(N)
% compute dual projection maps offline
    % use p to denote primal variables (such as rotation, translation, theta)
    t0 = tic;
    nrPrimalVars = 9 + 3 + N;
    p = mpvar('p',[nrPrimalVars,1]);
    r = p(1:9); R = reshape(r,[3,3]);
    r1 = R(:,1); r2 = R(:,2); r3 = R(:,3);
    t = p(10:12);
    theta = p(13:end);
    % build my constraints
    g_r = [1-r1'*r1;1-r2'*r2;1-r3'*r3;...
           r1'*r2;r2'*r3;r3'*r1;...
           hatmap(r1)*r2-r3;hatmap(r2)*r3-r1;hatmap(r3)*r1-r2];
    g_t = 1.0^2-t'*t;
    g_theta = [];
    for i = 1:N 
        g_theta =[g_theta; 1-theta(i)^2];
    end
    nrCons_r = length(g_r);
    nrCons_t = length(g_t);
    nrCons_theta = length(g_theta);

    % use d to denote dual variables (PSD matrices, lambdas)
    % f(p) - f_est - \sum lambda_ri*g_ri - \sum lambda_thetai*g_thetai = (basis_t' * Xt * basis_t) + (basis_p' * Xp * basis_p)
    % lambda_ri = d_ri' * basis_r
    % lambda_thetai = d_thetai' * basis_theta
    basis_r = monomials([t;theta],0:2);
    basis_theta = monomials([r;t],0:2);
    basis_t = monomials([r;theta],0:1);
    basis_p = [ monomials(p,0:1); kron(r,t); kron(theta,r); kron(theta,t) ];
    dim_Xt = length(basis_t);
    dim_Xp = length(basis_p);
    
    nrDualVars_r = length(basis_r) * nrCons_r;
    nrDualVars_theta = length(basis_theta) * nrCons_theta;
    nrDualVars_t = dim_Xt^2;
    nrDualVars_p = dim_Xp^2;
    nrDualVars_free = nrDualVars_r + nrDualVars_theta;
    nrDualVars_psd = nrDualVars_t + nrDualVars_p;
    nrDualVars = nrDualVars_free + nrDualVars_t + nrDualVars_p;
    % Build the affine subspace [A1, A2, A3] * d = b
    poly_fuse = [kron(g_r,basis_r);...
                 kron(g_theta,basis_theta);...
                 kron(g_t*basis_t,basis_t);...
                 kron(basis_p,basis_p)];

    assert(length(poly_fuse) == nrDualVars, 'Dimension is wrong.');

    A = poly_fuse.coefficient;
    A1 = A(:,1:nrDualVars_free);
    A2_raw = A(:,nrDualVars_free+1:nrDualVars_free+nrDualVars_t)';
    A3_raw = A(:,nrDualVars_free+nrDualVars_t+1:end)';
    nrTerms = size(A,1);
    nrTerms_full = nchoosek(nrPrimalVars + 2*2, 2*2);
    % convert to svec format to speed up
    dim_svect = dim_Xt * (dim_Xt + 1)/2;
    dim_svecp = dim_Xp * (dim_Xp + 1)/2;
    nrDualVars_svec = nrDualVars_free + dim_svect + dim_svecp;
    A2 = sparse(dim_svect,nrTerms);
    A3 = sparse(dim_svecp,nrTerms);
    for i = 1:size(A2_raw,2)
        A2(:,i) = svec(reshape(A2_raw(:,i),[dim_Xt,dim_Xt]));
    end
    for i = 1:size(A3_raw,2)
        A3(:,i) = svec(reshape(A3_raw(:,i),[dim_Xp,dim_Xp]));
    end
    A2 = A2'; A3 = A3';
    A = [A1, A2, A3];

    [inverse_A, inverse_preb] = get_dual_projection(A1,A2,A3);
    projection.type = 'mesh registration';
    projection.N = N;
    projection.A = A;
    projection.inv_A = inverse_A;
    projection.pre_b = inverse_preb;
    projection.degmat = poly_fuse.degmat;
    projection.varname = poly_fuse.varname;
    projection.nrDualVars_free = nrDualVars_free;
    projection.dim_svect = dim_svect;
    projection.dim_svecp = dim_svecp;
    projection.dim_Xt = dim_Xt;
    projection.dim_Xp = dim_Xp;

    t_total = toc(t0);

    fprintf('projection mesh registration: N=%d.\n',N);
    fprintf('# of Terms w/ reduction: \t %d\n',nrTerms);
    fprintf('# of Terms w/o reduction: \t %d\n',nrTerms_full);
    fprintf('# of Columns in A1: \t %d\n',size(A1,2));
    fprintf('# of Columns in A2: \t %d\n',size(A2,2));
    fprintf('# of Columns in A3: \t %d (Orthogonal)\n',size(A3,2));
    fprintf('Obtained projection map in %g seconds.\n',t_total);

    % save the projection map
    % assume this function is called from certifiable-perception/example_**
    % directories
    filename = sprintf('../solver/mesh_registration/projection_maps/projection_mr_%d.mat',N);
    save(filename,'projection');
end

function [proj_dual_A,proj_dual_pre_b] = get_dual_projection(A1,A2,A3)
    [nrTerms,~] = size(A1);
    A12 = [A1,A2];
    A = [A1,A2,A3];
    P = A3*A3';
    P_inv_diag = 1./diag(P);
    P_inv = spdiags(P_inv_diag,0,nrTerms,nrTerms);
    Z = speye(size(A12,2)) + A12'*P_inv*A12;
    M_inv = P_inv - P_inv * A12 * (Z \ A12') * P_inv;

    proj_dual_pre_b = A'*M_inv;
    proj_dual_A = speye(size(A,2)) - proj_dual_pre_b*A;
end