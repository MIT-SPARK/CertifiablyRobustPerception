function projection = get_dual_projection_shape_alignment_tf(N,varargin)
    % compute projection map offline, for shape alignment translation free

    params = inputParser;
    params.CaseSensitive = false;
    params.addParameter('s_ub',2.0,@(x) isscalar(x));
    params.parse(varargin{:});
    s_ub = params.Results.s_ub;

    fprintf('SATF dual projection: scale upperbound = %g.\n',s_ub);

    t0 = tic;
    nrPrimalVars = 6 + N;
    p = mpvar('p',[nrPrimalVars,1]);
    r = p(1:6); R = reshape(r,[2,3]);
    r1 = R(1,:)'; r2 = R(2,:)';
    theta = p(7:end);

    % build my constraints
    g_r_eq = [r1'*r1-r2'*r2;...
              r1'*r2];
    g_r_ineq = 2*s_ub^2 - r'*r;
    g_theta = [];
    for i = 1:N 
        g_theta =[g_theta; 1-theta(i)^2];
    end
    nrCons_r_eq = length(g_r_eq);
    nrCons_r_ineq = length(g_r_ineq);
    nrCons_theta = length(g_theta);

    % use d to denote dual variables (PSD matrices, lambdas)
    % f(p) - f_est - \sum lambda_ri*g_r_eqi - \sum lambda_thetai*g_thetai = (basis_r_ineq' * Xr * basis_r_ineq) + (basis_p' * Xp * basis_p)
    % lambda_ri = d_ri' * basis_r
    % lambda_thetai = d_thetai' * basis_theta
    basis_r_eq = monomials(theta,0:2);
    basis_theta = monomials(r,0:2);
    basis_r_ineq = monomials(theta,0:1);
    basis_p = [ monomials(p,0:1); kron(theta,r)];
    dim_Xr = length(basis_r_ineq);
    dim_Xp = length(basis_p);

    nrDualVars_r_eq = length(basis_r_eq) * nrCons_r_eq;
    nrDualVars_theta = length(basis_theta) * nrCons_theta;
    nrDualVars_r_ineq = dim_Xr^2;
    nrDualVars_p = dim_Xp^2;
    nrDualVars_free = nrDualVars_r_eq + nrDualVars_theta;
    nrDualVars_psd = nrDualVars_r_ineq + nrDualVars_p;
    nrDualVars = nrDualVars_free + nrDualVars_psd;
    % Build the affine subspace [A1, A2, A3] * d = b
    poly_fuse = [kron(g_r_eq,basis_r_eq);...
                 kron(g_theta,basis_theta);...
                 kron(g_r_ineq*basis_r_ineq,basis_r_ineq);...
                 kron(basis_p,basis_p)];

    assert(length(poly_fuse) == nrDualVars, 'Dimension is wrong.');

    A = poly_fuse.coefficient;
    A1 = A(:,1:nrDualVars_free);
    A2_raw = A(:,nrDualVars_free+1:nrDualVars_free+nrDualVars_r_ineq)';
    A3_raw = A(:,nrDualVars_free+nrDualVars_r_ineq+1:end)';
    nrTerms = size(A,1);
    nrTerms_full = nchoosek(nrPrimalVars + 2*2, 2*2);
    % convert to svec format to speed up
    dim_svecr = dim_Xr * (dim_Xr + 1)/2;
    dim_svecp = dim_Xp * (dim_Xp + 1)/2;
    nrDualVars_svec = nrDualVars_free + dim_svecr + dim_svecp;
    A2 = sparse(dim_svecr,nrTerms);
    A3 = sparse(dim_svecp,nrTerms);
    for i = 1:size(A2_raw,2)
        A2(:,i) = svec(reshape(A2_raw(:,i),[dim_Xr,dim_Xr]));
    end
    for i = 1:size(A3_raw,2)
        A3(:,i) = svec(reshape(A3_raw(:,i),[dim_Xp,dim_Xp]));
    end
    A2 = A2'; A3 = A3';
    A = [A1, A2, A3];

    [inverse_A, inverse_preb] = get_dual_projection(A1,A2,A3);
    projection.type = 'shape alignment translation free';
    projection.N = N;
    projection.A = A;
    projection.inv_A = inverse_A;
    projection.pre_b = inverse_preb;
    projection.degmat = poly_fuse.degmat;
    projection.varname = poly_fuse.varname;
    projection.nrDualVars_free = nrDualVars_free;
    projection.dim_svecr = dim_svecr;
    projection.dim_svecp = dim_svecp;
    projection.dim_Xr = dim_Xr;
    projection.dim_Xp = dim_Xp;

    t_total = toc(t0);

    fprintf('projection shape alignment translation free: N=%d.\n',N);
    fprintf('# of Terms w/ reduction: \t %d\n',nrTerms);
    fprintf('# of Terms w/o reduction: \t %d\n',nrTerms_full);
    fprintf('# of Columns in A1: \t %d\n',size(A1,2));
    fprintf('# of Columns in A2: \t %d\n',size(A2,2));
    fprintf('# of Columns in A3: \t %d (Orthogonal)\n',size(A3,2));
    fprintf('Obtained projection map in %g seconds.\n',t_total);

    % save the projection map
    % assume this function is called from certifiable-perception/example_**
    % directories
    filename = sprintf('../solver/shape_alignment_tf/projection_maps/projection_satf_%d.mat',N);
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