function solution = tls_point_cloud_registration_gnc(problem)
    % use GNC-TLS to solve point cloud registration
    % Heng Yang
    % 04/06/2020
    % copy data
    N = problem.N;
    cloudA = problem.cloudA;
    cloudB = problem.cloudB;
    noiseBoundSq = problem.noiseBoundSq;
    barc2 = 1; % since the residuals are normalized, barc2 is always 1
    translationBound = problem.translationBound;
    divFactor = 1.4;

    t0 = tic;
    allPoints = 1:N;  
    weights = ones(1,N);
    stopTh = 1e-16;
    maxSteps = 1e5;
    itr = 0;
    pre_f_cost = inf;
    f_cost = inf;
    cost_diff = inf;
    
    weightedCost=[];
    TLSCost = [];
    GNCCost=[];
    
    while itr < maxSteps && cost_diff > stopTh
        if max(abs(weights)) < 1e-6
            disp('Weights vanish, GNC failed.')
            break;
        end

        % fix weights and solve for s, t, and R using procrustes
        weighted_model = zeros(size(cloudA));
        weighted_scene = weighted_model;
        for i = 1:N
            weighted_model(:,i) = weights(i) * cloudA(:,i);
            weighted_scene(:,i) = weights(i) * cloudB(:,i);
        end
        weighted_model_center = sum(weighted_model,2) / sum(weights);
        weighted_scene_center = sum(weighted_scene,2) / sum(weights);
        
        demeaned_weighted_model = zeros(size(cloudA));
        demeaned_weighted_scene = demeaned_weighted_model;
        for i = 1:N
            demeaned_weighted_model(:,i) = (cloudA(:,i) - weighted_model_center) * weights(i)^(0.5);
            demeaned_weighted_scene(:,i) = (cloudB(:,i) - weighted_scene_center) * weights(i)^(0.5);
        end
        [f_cost, ~, transform] = procrustes(demeaned_weighted_scene', demeaned_weighted_model', ...
            'Scaling', false, ...
            'Reflection', false);
        R_est_cf = transform.T';
        s_est_cf = transform.b;
        t_est_cf = weighted_scene_center - s_est_cf * R_est_cf * weighted_model_center;

        % fix s, t, and R, solve for weights in closed form
        residuals = zeros(1,N);
        for i = 1:N
            residuals(i) = norm( cloudB(:,i) -  s_est_cf * R_est_cf * cloudA(:,i) - t_est_cf )^2 / noiseBoundSq;
        end
        
        if itr < 1
            maxResidual = max(residuals);
            mu = max(1 / ( 5 * maxResidual / barc2 - 1 ),1e-6); % make sure mu is postive
            % cprintf('Keywords', 'maxResidual=%g, set mu=%g.\n',maxResidual, mu);
        end
        
        
        th1 = (mu+1)/mu * barc2;
        th2 = (mu)/(mu+1) * barc2; % th1 > th2
        for i = 1:N
            if residuals(i) - th1 >= 0
                weights(i) = 0;
            elseif residuals(i) - th2 <= 0
                weights(i) = 1;
            else
                weights(i) = sqrt( barc2*mu*(mu+1)/residuals(i) ) - mu;
                assert(weights(i)>= 0 && weights(i) <=1, 'weights calculation wrong!');
            end
        end

        cost_diff = abs(f_cost - pre_f_cost);
        
        % increase mu
        mu = mu * divFactor;
        itr = itr + 1;
        pre_f_cost = f_cost;
    end
    
    R_est = R_est_cf;
    t_est = t_est_cf;
    s_est = s_est_cf;
    t_gnc = toc(t0);
    
    f_est = 0;
    for i = 1:N 
        residual = norm( cloudB(:,i) - R_est * cloudA(:,i) - t_est  )^2 / noiseBoundSq;
        f_est = f_est + min(barc2, residual);
    end

    solution.type = 'GNC-TLS';
    solution.weights = weights;
    theta_est = zeros(N,1);
    theta_est(weights > 0.5) = 1;
    theta_est(weights < 0.5) = -1;
    solution.theta_est = theta_est;
    solution.R_est = R_est;
    solution.t_est = t_est;
    solution.itr = itr;
    solution.divFactor = divFactor;
    solution.t_gnc = t_gnc;
    solution.f_est = f_est;
    solution.detectedOutliers = allPoints(weights<0.5);

    % print some info
    fprintf('============================== GNC-TLS ================================\n')
    fprintf('f_est = %g, divFactor=%g, itr=%d, t_gnc=%g[s].\n',f_est,divFactor,solution.itr,t_gnc);
    fprintf('=======================================================================\n')
    
end

