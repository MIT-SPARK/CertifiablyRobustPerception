function certification = tls_point_cloud_registration_certification_projection(problem,solution,varargin)
    % Given a solution to the point cloud registation problem
    % Certify the correctness of the solution
    % Compute a sub-optimality bound
    % Alternating projections to PSD Cone and Affine subspace
    % Use multipoly package to manipulate polynomials
    % Heng Yang
    % 04/07/2020

    params = inputParser;
    params.CaseSensitive = false;
    params.addParameter('method', 'Nesterov', ...
        @(x) ischar(x) && any(strcmpi({'Nesterov', 'Naive', 'D-R'}, x)));
    params.parse(varargin{:});
    method = params.Results.method;
    
    % copy problem data
    N = problem.N;
    cloudA = problem.cloudA;
    cloudB = problem.cloudB;
    noiseBoundSq = problem.noiseBoundSq;
    barc2 = 1; % since the residuals are normalized, barc2 is always 1
    translationBound = problem.translationBound;
    % copy solution data
    R_est = solution.R_est;
    r_est = R_est(:);
    t_est = solution.t_est;
    theta_est = ones(N,1);
    theta_est(solution.detectedOutliers) = -1.0;
    f_est = solution.f_est;

    %% compute projection map to the affine subspace
    projection = compute_dual_projection(problem,f_est);
    A = projection.A;
    b = projection.b;
    proj_dual_A = projection.inv_A;
    proj_dual_b = projection.inv_b;
    nrDualVars_free = projection.nrDualVars_free;
    dim_svect = projection.dim_svect;
    dim_svecp = projection.dim_svecp;
    dim_Xt = projection.dim_Xt;
    dim_Xp = projection.dim_Xp;
    
    %% start alternating projections
    alpha_psd = 2;
    alpha_affine = 1;
    maxIters = 1e4;
    xt_init = svec(speye(dim_Xt));
    xp_init = svec(speye(dim_Xp));
    xfree_init = zeros(nrDualVars_free,1);
    x_init = [xfree_init;xt_init;xp_init];
    x_affine = proj_dual_A * x_init + proj_dual_b;
    prev_x_affine = x_affine;
    z = x_affine;
    itr = 0;
    relDualityGapTol = 1e-2;
    suboptimality_traj = [];
    while itr < maxIters
        switch method
        case 'Nesterov'
            % compute suboptimality gap
            optInfo = compute_suboptimality(x_affine,nrDualVars_free,dim_svect,dim_svecp,f_est,N,translationBound);
            gamma = (itr-1)/(itr+2);
            v_affine = x_affine + gamma*(x_affine - prev_x_affine);
            prev_x_affine = x_affine;
            x_psd = dual_proj_psd(v_affine,nrDualVars_free,dim_svect,dim_svecp,alpha_psd);
            x_affine_1 = proj_dual_A * x_psd + proj_dual_b;
            x_affine = (1-alpha_affine) * x_affine + alpha_affine * x_affine_1;
        case 'Naive'
            optInfo = compute_suboptimality(x_affine,nrDualVars_free,dim_svect,dim_svecp,f_est,N,translationBound);
            x_psd = dual_proj_psd(x_affine,nrDualVars_free,dim_svect,dim_svecp,1);
            x_affine = proj_dual_A * x_psd + proj_dual_b;
        case 'D-R'
            % gamma = (itr-1)/(itr+2);
            gamma = 1;
            x_psd = dual_proj_psd(z,nrDualVars_free,dim_svect,dim_svecp,1);
            x_affine = proj_dual_A * (2*x_psd - z) + proj_dual_b;
            z = z + 1 * (x_affine - x_psd);
            optInfo = compute_suboptimality(x_affine,nrDualVars_free,dim_svect,dim_svecp,f_est,N,translationBound);
        otherwise 
            error('Unsupported projection method.');    
        end

        suboptimality_traj = [suboptimality_traj,optInfo.relDualityGap];
        if optInfo.relDualityGap < relDualityGapTol
            fprintf('Relative duality gap below %g percent at %d iteration.\n',relDualityGapTol*100,itr);
            break;
        end
        if rem(itr,1000) == 0
            fprintf('Itr: %d, relDualityGap = %g, minEigs = %g, %g.\n',itr,optInfo.relDualityGap,optInfo.min_eig_Xt,optInfo.min_eig_Xp);
            % disp(norm(A*x_affine - b))
        end
        itr = itr + 1;
    end
    figure;
    plot(suboptimality_traj,'linewidth',2);
    set(gca,'yscale','log')

    certification.projection = projection;
    certification.optInfo = optInfo;
    
    
end

function optInfo = compute_suboptimality(x,d1,d2,d3,f_est,N,translationBound)
    Xt = smat(x(d1+1:d1+d2));
    Xp = smat(x(d1+d2+1:end));
    min_eig_Xt = compute_min_eig(Xt);
    min_eig_Xp = compute_min_eig(Xp);
    absDualityGap = (4+N)*min_eig_Xt + (4+4*translationBound^2 + (4+translationBound^2) * N) * min_eig_Xp;
    absDualityGap = abs( absDualityGap );
    relDualityGap = absDualityGap / abs(f_est);
    % fprintf('min eigs: %g, %g; dualityGap: %g, %g.\n',min_eig_Xt,min_eig_Xp,absDualityGap,relDualityGap);
    optInfo.min_eig_Xt = min_eig_Xt;
    optInfo.min_eig_Xp = min_eig_Xp;
    optInfo.absDualityGap = absDualityGap;
    optInfo.relDualityGap = relDualityGap;
end

function min_eig = compute_min_eig(A)
    A = (A + A')/2;
    min_eig = min(eig(A));
end

function x_proj = dual_proj_psd(x,d1,d2,d3,alpha_psd)
    x = full(x);
    x_free = x(1:d1);
    x_t = x(d1+1:d1+d2);
    x_p = x(d1+d2+1:end);
    x_t_psd = svec(simpleNearestSPD(smat(x_t)));
    x_p_psd = svec(simpleNearestSPD(smat(x_p)));
    x_proj = [x_free;x_t_psd;x_p_psd];
    x_proj = (1-alpha_psd) * x + alpha_psd * x_proj;
end

function projection = compute_dual_projection(problem,f_est)
    t00 = tic;
    N = problem.N;
    noiseBoundSq = problem.noiseBoundSq;
    cloudA = problem.cloudA;
    cloudB = problem.cloudB;
    barc2 = 1.0;
    
    % load or compute projection matrix
    filename = sprintf('../solver/point_cloud_registration/projection_maps/projection_pcr_%d.mat',N);
    if isfile(filename)
        t0 = tic; load(filename); t_load = toc(t0);
        fprintf('Load projection file in %g seconds.\n',t_load);
    else
        disp('No projection file exists, computing now...')
        projection = get_dual_projection_point_cloud_registration(N);
    end
    % use p to denote primal variables (such as rotation, translation, theta)
    nrPrimalVars = 9 + 3 + N;
    p = mpvar('p',[nrPrimalVars,1]);
    r = p(1:9); R = reshape(r,[3,3]);
    t = p(10:12);
    theta = p(13:end);
    % build my objective function
    residuals = {};
    for i = 1:N
        ai = cloudA(:,i);
        bi = cloudB(:,i);
        residual = bi'*bi + ai'*ai + t'*t - 2*bi'*R*ai - 2*bi'*t + 2*t'*R*ai;
        residuals{end+1} = residual / noiseBoundSq;
    end
    f_cost = 0;
    for i = 1:N
        fi = (1 + theta(i))/2 * residuals{i} + (1 - theta(i))/2 * barc2;
        f_cost = f_cost + fi;
    end

    g_cost = f_cost - f_est;
    nrTerms = size(projection.degmat,1);
    b = sparse(nrTerms,1);
    [~,IA,IB] = intersect(projection.degmat,g_cost.degmat,'rows','stable');
    b(IA) = g_cost.coefficient(IB);
    projection.b = b;
    projection.inv_b = projection.pre_b * b;
    t_proj = toc(t00);
    fprintf('Compute projection map for this problem in %g seconds.\n',t_proj);
end





% % Build the complementary slackness condition 
% basis_t_est = subs(basis_t,[r;theta],[r_est;theta_est]);
% basis_t_est = basis_t_est.coefficient';
% basis_p_est = subs(basis_p,p,[r_est;t_est;theta_est]);
% basis_p_est = basis_p_est.coefficient';
% At = sparse(dim_svect,dim_Xt);
% Ap = sparse(dim_svecp,dim_Xp);
% for i = 1:dim_Xt
%     temp = kron(sparse(1,i,1,1,dim_Xt), basis_t_est);
%     temp = (temp + temp')/2;
%     At(:,i) = svec(temp); 
% end
% for i = 1:dim_Xp
%     temp = kron(sparse(1,i,1,1,dim_Xp), basis_p_est);
%     temp = (temp + temp')/2;
%     Ap(:,i) = svec(temp); 
% end
% At = At'; Ap = Ap';