function certification = tls_shape_alignment_tf_certification_hybrid(problem,solution,varargin)
    % Certify shape alignment translation free
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

    global N f_est s_ub;
    N = problem.N;
    f_est = solution.f_est;
    s_ub = problem.scaleBound(2);

    fprintf('SATF certification: scale upperbound=%g.\n',s_ub);

    relDualityGapTol = 1e-2;
    certification.relDualityGapTol = relDualityGapTol;
    %% first find a chordal sparse SOS certificate
    chordal_certificate = find_chordal_sos_certificate_shape_alignment_tf(problem,solution);

    certification.chordal_certificate = chordal_certificate;

    if (chordal_certificate.relDualityGap > relDualityGapTol) || fixIters
        fprintf('======================= Douglas-Rachford Splitting ======================\n')
        t0 = tic;
        fprintf('Run first-order method: maxIters=%d, fixIters=%s, plotSuboptTraj=%s.\n',maxIters,logical2str(fixIters),logical2str(plotSuboptTraj));
        %% now use the solution from chordal SOS as an initial guess
        global nrDualVars_free dim_svecr dim_svecp;
        nrDualVars_free = chordal_certificate.nrDualVars_free;
        dim_svecr = chordal_certificate.dim_svecr;
        dim_svecp = chordal_certificate.dim_svecp;
        x0 = chordal_certificate.x0;

        % get the projection map
        projection = compute_dual_projection(problem,f_est);
        global proj_dual_A proj_dual_b
        proj_dual_A = projection.inv_A;
        proj_dual_b = projection.inv_b;

        %% start alternating projections from the initial guess
        z = x0;
        itr = 0;
        suboptimality_traj = [chordal_certificate.relDualityGap];
        gamma = 2;
        while itr < maxIters
            x_psd = dual_proj_psd(z);
            x_affine = proj_dual_A * (2*x_psd - z) + proj_dual_b;
            z = z + gamma * (x_affine - x_psd);
            optInfo = compute_suboptimality(x_affine);

            suboptimality_traj = [suboptimality_traj,optInfo.relDualityGap];

            if (optInfo.relDualityGap < relDualityGapTol) && (~fixIters)
                fprintf('Relative duality gap below %g percent at %d iteration.\n',relDualityGapTol*100,itr);
                break;
            end
            if rem(itr,1000) == 0
                fprintf('Itr: %d, relDualityGap = %g, minEigs = %g, %g.\n',...
                    itr,optInfo.relDualityGap,optInfo.min_eig_Xr,optInfo.min_eig_Xp);
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

        [bestRelDualityGap,idx] = min(suboptimality_traj);

        if bestRelDualityGap < relDualityGapTol
            fprintf('Global Optimality Certified in %g[s], bestRelDualityGap=%g.\n',...
                time_projection,bestRelDualityGap);
            certified = true;
        else
            fprintf('Cannot certify global optimality in %d iterations (%g[s]). Best suboptimality: %g.\n',...
                maxIters,time_projection,bestRelDualityGap);
            certified = false;
        end
        certification.nrDRSIters = itr;
        certification.proj_opt_info = optInfo;
        certification.projection = projection;
        certification.time_projection = time_projection;
        certification.time_certification = time_projection + chordal_certificate.time_chordal_sdp;
        certification.certified = certified;
        certification.suboptimality_traj = suboptimality_traj;
        certification.bestRelDualityGap = bestRelDualityGap;
    else
        disp('Chordal sos certificate is sufficient to guarantee global optimality.')
        certification.nrDRSIters = 0;
        certification.time_certification = chordal_certificate.time_chordal_sdp;
        certification.bestRelDualityGap = chordal_certificate.relDualityGap;
    end
    fprintf('=========================================================================\n')
end

function optInfo = compute_suboptimality(x)
    global nrDualVars_free dim_svecr dim_svecp f_est s_ub N;
    Xr = smat(x(nrDualVars_free+1:nrDualVars_free+dim_svecr));
    Xp = smat(x(nrDualVars_free+dim_svecr+1:end));
    min_eig_Xr = compute_min_eig(Xr);
    min_eig_Xp = compute_min_eig(Xp);
    absDualityGap = (1+N)*min_eig_Xr + (1+2*s_ub^2)*(1+N) * min_eig_Xp;
    absDualityGap = abs( absDualityGap );
    relDualityGap = absDualityGap / abs(f_est);
    optInfo.min_eig_Xr = min_eig_Xr;
    optInfo.min_eig_Xp = min_eig_Xp;
    optInfo.absDualityGap = absDualityGap;
    optInfo.relDualityGap = relDualityGap;
end

function min_eig = compute_min_eig(A)
    A = (A + A')/2;
    min_eig = min(eig(A));
end

function x_proj = dual_proj_psd(x)
    global nrDualVars_free dim_svecr dim_svecp;
    x = full(x);
    x_free = x(1:nrDualVars_free);
    x_r = x(nrDualVars_free+1:nrDualVars_free+dim_svecr);
    x_p = x(nrDualVars_free+dim_svecr+1:end);
    x_r_psd = svec(simpleNearestSPD(smat(x_r)));
    x_p_psd = svec(simpleNearestSPD(smat(x_p)));
    x_proj = [x_free;x_r_psd;x_p_psd];
end

function projection = compute_dual_projection(problem,f_est)
    t00 = tic;
    N = problem.N;
    z = problem.z;
    B = problem.B;
    noiseBoundSq = problem.noiseBoundSq;
    barc2 = 1; % since the residuals are normalized, barc2 is always 1
    
    % load or compute projection matrix
    filename = sprintf('../solver/shape_alignment_tf/projection_maps/projection_satf_%d.mat',N);
    if isfile(filename)
        t0 = tic; load(filename); t_load = toc(t0);
        fprintf('Load projection file in %g seconds.\n',t_load);
    else
        disp('No projection file exists, computing now...')
        projection = get_dual_projection_shape_alignment_tf(N);
    end
    % use p to denote primal variables (such as rotation, translation, theta)
    nrPrimalVars = 6 + N;
    p = mpvar('p',[nrPrimalVars,1]);
    r = p(1:6); R = reshape(r,[2,3]);
    theta = p(7:end);
    % build my objective function
    residuals = {};
    for i = 1:N
        zi = z(:,i);
        Bi = B(:,i);
        resVec = zi - R*Bi;
        
        residuals{end+1} = resVec'*resVec / noiseBoundSq;
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