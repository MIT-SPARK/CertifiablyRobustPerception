function certification = tls_single_rotation_averaging_certification_hybrid(problem,solution,varargin)
    % certify single rotation averaging
    % Heng Yang
    % 04/20/2020

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

    global N f_est vp_est;
    N = problem.N;
    f_est = solution.f_est;
    R_est = solution.R_est; r_est = R_est(:);
    theta_est = solution.theta_est;
    relDualityGapTol = 1e-2;
    PDSlacknessTol = 1e-2;
    certification.relDualityGapTol = relDualityGapTol;
    certification.PDSlacknessTol = PDSlacknessTol;
    %% first find a chordal sparse SOS certificate
    chordal_certificate = find_chordal_sos_certificate_single_rotation_averaging(problem,solution);
    vp_est = chordal_certificate.vp_est;
    certification.chordal_certificate = chordal_certificate;
    if (chordal_certificate.relDualityGap > relDualityGapTol) || fixIters
        fprintf('======================= Douglas-Rachford Splitting ======================\n')
        fprintf('Run first-order method: maxIters=%d, fixIters=%s, plotSuboptTraj=%s.\n',maxIters,logical2str(fixIters),logical2str(plotSuboptTraj));
        t0 = tic;
        global nrDualVars_free dim_svecp;
        nrDualVars_free = chordal_certificate.nrDualVars_free;
        dim_svecp = chordal_certificate.dim_svecp;
        x0 = chordal_certificate.x0;

        % get the projection map to affine subspace
        projection = compute_dual_projection(problem,f_est);
        global proj_dual_A proj_dual_b;
        proj_dual_A = projection.inv_A;
        proj_dual_b = projection.inv_b;

        % start alternating D-R 
        z = x0;
        itr = 0;
        gamma = 2;
        suboptimality_traj = [chordal_certificate.relDualityGap];
        PDSlackness_traj = [chordal_certificate.PDSlackness];
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
                fprintf('Itr: %d, relDualityGap = %g, minEigs = %g, PDSlackness = %g.\n',itr,optInfo.relDualityGap,optInfo.min_eig_Xp,optInfo.PDSlackness);
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

        time_projection = toc(t0);
        [bestRelDualityGap, idx] = min(suboptimality_traj);
        bestPDSlackness = PDSlackness_traj(idx);
        if optInfo.relDualityGap < relDualityGapTol
            fprintf('Global Optimality Certified in %g[s], relDualityGap = %g, PDSlackness = %g.\n',time_projection,bestRelDualityGap,bestPDSlackness);
            certified = true;
        else
            fprintf('Cannot certify global optimality in %d iterations (%g[s]). Best suboptimality: %g, best PDSlackness: %g.\n',maxIters,time_projection,bestRelDualityGap,bestPDSlackness);
            certified = false;
        end
        certification.proj_opt_info = optInfo;
        certification.projection = projection;
        certification.time_projection = time_projection;
        certification.time_certification = time_projection + chordal_certificate.time_chordal_sdp;
        certification.certified = certified;
        certification.suboptimality_traj = suboptimality_traj;
        certification.PDSlackness_traj = PDSlackness_traj;
        certification.bestRelDualityGap = bestRelDualityGap;
        certification.bestPDSlackness = bestPDSlackness;
    else
        disp('Chordal sos certificate is sufficient to guarantee global optimality.')
        certification.time_certification = chordal_certificate.time_chordal_sdp;
    end
    fprintf('=========================================================================\n')
end

function optInfo = compute_suboptimality(x)
    global nrDualVars_free dim_svecp f_est N vp_est;
    Xp = smat(x(nrDualVars_free+1:end));
    min_eig_Xp = compute_min_eig(Xp);
    absDualityGap = (4*N+4) * min_eig_Xp;
    absDualityGap = abs( absDualityGap );
    relDualityGap = absDualityGap / abs(f_est);
    % fprintf('min eigs: %g, %g; dualityGap: %g, %g.\n',min_eig_Xt,min_eig_Xp,absDualityGap,relDualityGap);
    optInfo.min_eig_Xp = min_eig_Xp;
    optInfo.absDualityGap = absDualityGap;
    optInfo.relDualityGap = relDualityGap;
    PDSlackness = norm(Xp * vp_est);
    optInfo.PDSlackness = PDSlackness;
end

function min_eig = compute_min_eig(A)
    A = (A + A')/2;
    min_eig = min(eig(A));
end

function x_proj = dual_proj_psd(x)
    global nrDualVars_free dim_svecp
    x = full(x);
    x_free = x(1:nrDualVars_free);
    x_p = x(nrDualVars_free+1:end);
    x_p_psd = svec(simpleNearestSPD(smat(x_p)));
    x_proj = [x_free;x_p_psd];
    % x_proj = (1-alpha_psd) * x + alpha_psd * x_proj;
end

function projection = compute_dual_projection(problem,f_est)
    t00 = tic;
    N = problem.N;
    noiseBoundSq = problem.noiseBoundSq;
    R_measurements = problem.R_measurements;
    barc2 = 1.0;
    
    % load or compute projection matrix
    filename = sprintf('../solver/single_rotation_averaging/projection_maps/projection_sra_%d.mat',N);
    if isfile(filename)
        t0 = tic; load(filename); t_load = toc(t0);
        fprintf('Load projection file in %g seconds.\n',t_load);
    else
        disp('No projection file exists, computing now...')
        projection = get_dual_projection_single_rotation_averaging(N);
    end
    % use p to denote primal variables (such as rotation, translation, theta)
    nrPrimalVars = 9 + N;
    p = mpvar('p',[nrPrimalVars,1]);
    r = p(1:9); R = reshape(r,[3,3]);
    theta = p(10:end);
    % build my objective function
    residuals = {};
    for i = 1:N
        ri = reshape(R_measurements(:,:,i),[9,1]);
        residual = (r-ri)'*(r-ri);
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


