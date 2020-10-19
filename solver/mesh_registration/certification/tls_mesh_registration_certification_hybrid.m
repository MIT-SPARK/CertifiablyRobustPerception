function certification = tls_mesh_registration_certification_hybrid(problem,solution,varargin)
    % certify the correctness of mesh registration solution 

    params = inputParser;
    params.CaseSensitive = false;
    params.addParameter('plotSuboptTraj',true,...
        @(x) islogical(x));
    params.addParameter('maxIters',2e3,@(x) isscalar(x));
    params.addParameter('fixIters',false,@(x) islogical(x));
    params.parse(varargin{:});
    plotSuboptTraj = params.Results.plotSuboptTraj;
    maxIters = params.Results.maxIters;
    fixIters = params.Results.fixIters; % fix number of iterations, run until maxIters

    global N translationBound f_est vt_est vp_est;
    N = problem.N;
    translationBound = problem.translationBound;
    f_est = solution.f_est;

    relDualityGapTol = 1e-2;
    PDSlacknessTol = 1e-2;
    certification.relDualityGapTol = relDualityGapTol;
    certification.PDSlacknessTol = PDSlacknessTol;
    %% first find a chordal sparse SOS certificate
    chordal_certificate = find_chordal_sos_certificate_mesh_registration(problem,solution);
    vt_est = chordal_certificate.vt_est;
    vp_est = chordal_certificate.vp_est;

    certification.chordal_certificate = chordal_certificate;
    if (chordal_certificate.relDualityGap > relDualityGapTol) || fixIters
        fprintf('======================= Douglas-Rachford Splitting ======================\n')
        t0 = tic;
        fprintf('Run first-order method: maxIters=%d, fixIters=%s, plotSuboptTraj=%s.\n',maxIters,logical2str(fixIters),logical2str(plotSuboptTraj));
        %% now use the solution from chordal SOS as an initial guess
        global nrDualVars_free dim_svect dim_svecp;
        nrDualVars_free = chordal_certificate.nrDualVars_free;
        dim_svect = chordal_certificate.dim_svect;
        dim_svecp = chordal_certificate.dim_svecp;
        x0 = chordal_certificate.x0;

        % get the projection map
        projection = compute_dual_projection(problem,f_est);
        A = projection.A;
        b = projection.b;
        global proj_dual_A proj_dual_b
        proj_dual_A = projection.inv_A;
        proj_dual_b = projection.inv_b;

        %% start alternating projections from the initial guess
        z = x0;
        itr = 0;
        suboptimality_traj = [chordal_certificate.relDualityGap];
        PDSlackness_traj = [chordal_certificate.PDSlackness];
        gamma = 2;
        while itr < maxIters
            x_psd = dual_proj_psd(z);
            x_affine = proj_dual_A * (2*x_psd - z) + proj_dual_b;
            z = z + gamma * (x_affine - x_psd);
            optInfo = compute_suboptimality(x_affine);

            suboptimality_traj = [suboptimality_traj,optInfo.relDualityGap];
            PDSlackness_traj = [PDSlackness_traj,optInfo.PDSlackness];

            if (optInfo.relDualityGap < relDualityGapTol) && (~fixIters)
                fprintf('Relative duality gap below %g percent at %d iteration.\n',relDualityGapTol*100,itr);
                break;
            end
            if rem(itr,1000) == 0
                fprintf('Itr: %d, relDualityGap = %g, minEigs = %g, %g, PDSlackness = %g, %g.\n',...
                    itr,optInfo.relDualityGap,optInfo.min_eig_Xt,optInfo.min_eig_Xp,optInfo.PDSlackness(1),optInfo.PDSlackness(2));
            end
            itr = itr + 1;
        end
        if plotSuboptTraj
            figure;
            plot(suboptimality_traj,'linewidth',2);
            set(gca,'yscale','log');
            xlabel('Iteration');
            ylabel('Suboptimality');
            grid on
        end

        [bestRelDualityGap,idx] = min(suboptimality_traj);
        bestPDSlackness = PDSlackness_traj(:,idx);
        time_projection = toc(t0);
        if optInfo.relDualityGap < relDualityGapTol
            fprintf('Global Optimality Certified in %g[s], bestRelDualityGap=%g, bestPDSlackness=%g, %g.\n',...
                time_projection,bestRelDualityGap,bestPDSlackness(1),bestPDSlackness(2));
            certified = true;
        else
            fprintf('Cannot certify global optimality in %d iterations (%g[s]). Best suboptimality: %g, Best PDSlackness: %g, %g.\n',...
                maxIters,time_projection,bestRelDualityGap,bestPDSlackness(1),bestPDSlackness(2));
            certified = false;
        end
        certification.proj_opt_info = optInfo;
        certification.projection = projection;
        certification.time_projection = time_projection;
        certification.time_certification = time_projection + chordal_certificate.time_chordal_sdp;
        certification.certified = certified;
        certification.suboptimality_traj = suboptimality_traj;
%         certification.PDSlackness_traj = PDSlackness_traj;
        certification.bestRelDualityGap = bestRelDualityGap;
        certification.bestPDSlackness = bestPDSlackness;
    else
        disp('Chordal sos certificate is sufficient to guarantee global optimality.')
        certification.time_certification = chordal_certificate.time_chordal_sdp;
    end
    fprintf('=========================================================================\n')
end

function optInfo = compute_suboptimality(x)
    global nrDualVars_free dim_svect dim_svecp translationBound N f_est vt_est vp_est;
    Xt = smat(x(nrDualVars_free+1:nrDualVars_free+dim_svect));
    Xp = smat(x(nrDualVars_free+dim_svect+1:end));
    min_eig_Xt = compute_min_eig(Xt);
    min_eig_Xp = compute_min_eig(Xp);
    absDualityGap = (4+N)*min_eig_Xt + (4+4*translationBound^2 + (4+translationBound^2) * N) * min_eig_Xp;
    absDualityGap = abs( absDualityGap );
    relDualityGap = absDualityGap / abs(f_est);
    optInfo.min_eig_Xt = min_eig_Xt;
    optInfo.min_eig_Xp = min_eig_Xp;
    optInfo.absDualityGap = absDualityGap;
    optInfo.relDualityGap = relDualityGap;
    PDSlackness = [norm(Xt*vt_est); norm(Xp*vp_est)];
    optInfo.PDSlackness = PDSlackness;
end

function min_eig = compute_min_eig(A)
    A = (A + A')/2;
    min_eig = min(eig(A));
end

function x_proj = dual_proj_psd(x)
    global nrDualVars_free dim_svect dim_svecp;
    x = full(x);
    x_free = x(1:nrDualVars_free);
    x_t = x(nrDualVars_free+1:nrDualVars_free+dim_svect);
    x_p = x(nrDualVars_free+dim_svect+1:end);
    x_t_psd = svec(simpleNearestSPD(smat(x_t)));
    x_p_psd = svec(simpleNearestSPD(smat(x_p)));
    x_proj = [x_free;x_t_psd;x_p_psd];
end

function projection = compute_dual_projection(problem,f_est)
    t00 = tic;
    N = problem.N;
    normalM = problem.normalM;
    normalP = problem.normalP;
    pointM = problem.pointM;
    pointP = problem.pointP;
    pointNoiseBoundSq = problem.pointNoiseBoundSq;
    normalNoiseBoundSq = problem.normalNoiseBoundSq;
    barc2 = 1; % since the residuals are normalized, barc2 is always 1
    translationBound = problem.translationBound;
    
    % load or compute projection matrix
    filename = sprintf('../solver/mesh_registration/projection_maps/projection_mr_%d.mat',N);
    if isfile(filename)
        t0 = tic; load(filename); t_load = toc(t0);
        fprintf('Load projection file in %g seconds.\n',t_load);
    else
        disp('No projection file exists, computing now...')
        projection = get_dual_projection_mesh_registration(N);
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
        pointM_i = pointM(:,i);
        pointP_i = pointP(:,i);
        normalM_i = normalM(:,i);
        normalP_i = normalP(:,i);
        % note the objective function is different
        residual_point_i = ( normalM_i' * (R*pointP_i - pointM_i - t) )^2;
        residual_normal_i = normalP_i'*normalP_i + normalM_i'*normalM_i - 2*normalP_i'*R'*normalM_i;
        
        residuals{end+1} = (residual_point_i/pointNoiseBoundSq + residual_normal_i/normalNoiseBoundSq)/2;
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