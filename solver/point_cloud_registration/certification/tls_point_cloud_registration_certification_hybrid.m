function certification = tls_point_cloud_registration_certification_hybrid(problem,solution,varargin)
    % Given a solution to the point cloud registation problem
    % Certify the correctness of the solution
    % Compute a sub-optimality bound
    % Alternating projections to PSD Cone and Affine subspace
    % Use multipoly package to manipulate polynomials
    % Heng Yang
    % 04/12/2020

    params = inputParser;
    params.CaseSensitive = false;
    params.addParameter('method', 'D-R', ...
        @(x) ischar(x) && any(strcmpi({'Nesterov', 'Naive', 'D-R', 'PDA', 'PDAA'}, x)));
    params.addParameter('plotSuboptTraj',true,...
        @(x) islogical(x));
    params.addParameter('maxIters',2e3,@(x) isscalar(x));
    params.addParameter('fixIters',false,@(x) islogical(x));
    params.parse(varargin{:});
    method = params.Results.method;
    plotSuboptTraj = params.Results.plotSuboptTraj;
    maxIters = params.Results.maxIters;
    fixIters = params.Results.fixIters; % fix number of iterations, run until maxIters

    global N translationBound f_est;
    N = problem.N;
    translationBound = problem.translationBound;
    f_est = solution.f_est;
    relDualityGapTol = 1e-2;
    %% first find a chordal sparse SOS certificate
    chordal_certificate = find_chordal_sos_certificate_point_cloud_registration(problem,solution);
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
        alpha_psd = 2;
        alpha_affine = 1;
        x_affine = proj_dual_A * x0 + proj_dual_b;
        prev_x_affine = x_affine;
        itr = 0;
        suboptimality_traj = [chordal_certificate.relDualityGap];

        % for PDA
        xbar = x0;
        y = zeros(length(xbar),1);
        x = x0;
        prev_x = x0;
        pda_tau = 1;
        pda_sigma = 1;

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
            case 'D-R' % Douglas-Rachford splitting
                % gamma = (itr-1)/(itr+2);
                gamma = 2;
                x_psd = dual_proj_psd(z,nrDualVars_free,dim_svect,dim_svecp,1);
                x_affine = proj_dual_A * (2*x_psd - z) + proj_dual_b;
                z = z + gamma * (x_affine - x_psd);
                optInfo = compute_suboptimality(x_affine,nrDualVars_free,dim_svect,dim_svecp,f_est,N,translationBound);
            case 'PDA' % primal-dual algorithm (Chambolle and Pock 2011, A First-Order Primal-Dual Algorithm for Convex Problems with Applications to Imaging)
                pda_theta = 1;
                [y,x_affine] = pda_proj_affine_conjugate(xbar,y,pda_sigma);
                prev_x = x;
                x = pda_proj_psd(prev_x,y,pda_tau);
                xbar = x + pda_theta * (x - prev_x);
                optInfo = compute_suboptimality(x_affine,nrDualVars_free,dim_svect,dim_svecp,f_est,N,translationBound);    
            case 'PDAA' % primal-dual algorithm (Chambolle and Pock 2011, A First-Order Primal-Dual Algorithm for Convex Problems with Applications to Imaging)
                gamma = 0.001;
                [y,x_affine] = pda_proj_affine_conjugate(xbar,y,pda_sigma);
                prev_x = x;
                x = pda_proj_psd(prev_x,y,pda_tau);
                % update parameters
                pda_theta = 1 / sqrt(1+2*gamma*pda_tau);
                pda_tau = pda_theta * pda_tau;
                pda_sigma = pda_sigma / pda_theta;
                xbar = x + pda_theta * (x - prev_x);
                optInfo = compute_suboptimality(x_affine,nrDualVars_free,dim_svect,dim_svecp,f_est,N,translationBound);  
            otherwise 
                error('Unsupported projection method.');    
            end

            suboptimality_traj = [suboptimality_traj,optInfo.relDualityGap];

            if (optInfo.relDualityGap < relDualityGapTol) && (~fixIters)
                fprintf('Relative duality gap below %g percent at %d iteration.\n',relDualityGapTol*100,itr);
                break;
            end

            if rem(itr,1000) == 0
                fprintf('Itr: %d, relDualityGap = %g, minEigs = %g, %g.\n',itr,optInfo.relDualityGap,optInfo.min_eig_Xt,optInfo.min_eig_Xp);
                % disp(norm(A*x_affine - b))
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
        if optInfo.relDualityGap < relDualityGapTol
            fprintf('Global Optimality Certified in %g[s].\n',time_projection);
            certified = true;
        else
            fprintf('Cannot certify global optimality in %d iterations (%g[s]). Best suboptimality: %g percent.\n',maxIters,time_projection,min(suboptimality_traj)*100);
            certified = false;
        end
        certification.proj_opt_info = optInfo;
        certification.projection = projection;
        certification.time_projection = time_projection;
        certification.time_certification = time_projection + chordal_certificate.time_chordal_sdp;
        certification.certified = certified;
        certification.suboptimality_traj = suboptimality_traj;
        certification.bestRelDualityGap = min(suboptimality_traj);
    else
        disp('Chordal sos certificate is sufficient to guarantee global optimality.')
        certification.time_certification = chordal_certificate.time_chordal_sdp;
    end

    fprintf('=========================================================================\n')

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

function [y_next,x_affine] = pda_proj_affine_conjugate(xbar,y,sigma)
    global proj_dual_A proj_dual_b
    temp = y + sigma * xbar;
    conjugate = proj_dual_A * temp + sigma * proj_dual_b;
    y_next = temp - conjugate;
    x_affine = conjugate / sigma;
end

function x_next = pda_proj_psd(x,y,tau)
    temp = x - tau * y;
    global nrDualVars_free dim_svect dim_svecp
    x_next = dual_proj_psd(temp,nrDualVars_free,dim_svect,dim_svecp,1);
end
