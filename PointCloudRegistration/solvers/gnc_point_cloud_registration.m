function solution = gnc_point_cloud_registration(problem)
%% use graduated non-convexity to solve robust point cloud registration
%% Heng Yang
%% June 25, 2021

    N               = problem.N;
    cloudA          = problem.cloudA;
    cloudB          = problem.cloudB;
    noiseBoundSq    = problem.noiseBoundSq;
    barc2           = 1; % since the residuals are normalized, barc2 is always 1
    

    t0              = tic;
    fprintf('\n\n===================================== GNC-TLS =======================================\n')

    allPoints       = 1:N;  
    weights         = ones(1,N);
    stopTh          = 1e-20;
    divFactor       = 1.4;
    maxSteps        = 1e5;
    itr             = 0;
    pre_f_cost      = 1e6;
    f_cost          = 1e6;
    cost_diff       = 1e6;

    fprintf('epsilon: %3.2e, div: %1.1f, maxiters: %3.0e.\n',stopTh,divFactor,maxSteps);
    fprintf('-------------------------------------------------------------------------------------\n')
    fprintf(' itr |    obj        delobj  |    mu    |   sumw   |  sumout  |   maxres   |   gap   |\n')
    
    breakyes        = 0;
    while (true) 
        if max(abs(weights)) < 1e-6
            msg      = 'GNC encounters numerical issues, the solution is likely to be wrong.';
            breakyes = 3;
        end
        if itr == maxSteps
            breakyes = 2;
            msg      = 'Maximum iterations reached.';
        end
        if cost_diff < stopTh
            breakyes = 1;
            msg      = sprintf('GNC converged %3.2e < %3.2e.',cost_diff,stopTh);
        end
        if breakyes > 0
            fprintf('%s\n',msg);
            break
        end

        % fix weights and solve for s, t, and R using procrustes
        wsqrt    = sqrt(weights);
        A_center = (cloudA * weights') / sum(weights);
        B_center = (cloudB * weights') / sum(weights);
        
        demeaned_weighted_A = wsqrt .* (cloudA - A_center);
        demeaned_weighted_B = wsqrt .* (cloudB - B_center);

        [f_cost, ~, transform] = procrustes(demeaned_weighted_B', demeaned_weighted_A', ...
            'Scaling', false, ...
            'Reflection', false);
        R_est_cf = transform.T';
        s_est_cf = transform.b;
        t_est_cf = B_center - s_est_cf * R_est_cf * A_center;

        % fix s, t, and R, solve for weights in closed form
        clouddiff = cloudB - s_est_cf * R_est_cf * cloudA - t_est_cf;
        residuals = sum(clouddiff .^ 2,1) / noiseBoundSq;
        maxResidual        = max(residuals);
        f_cost_TLS         = sum(min(residuals,barc2));
        if itr < 1
            mu = max(1 / ( 5 * maxResidual / barc2 - 1 ),1e-6); % make sure mu is postive
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

        fprintf('%4d | %3.4e %3.4e | %3.2e | %3.2e | %3.2e | %3.4e | %1.1e |\n',...
            itr,f_cost_TLS,cost_diff,mu,sum(weights),N-sum(weights),maxResidual,0);
        
        % increase mu
        mu          = mu * divFactor;
        itr         = itr + 1;
        pre_f_cost  = f_cost;
    end
    
    R_est = R_est_cf;
    t_est = t_est_cf;
    s_est = s_est_cf;
    time_gnc = toc(t0);

    clouddiff = cloudB - s_est * R_est * cloudA - t_est;
    residuals = sum(clouddiff .^ 2,1) / noiseBoundSq;
    
    f_est     = sum(min(residuals,barc2));

    solution.type           = 'GNC-TLS';
    solution.weights        = weights;
    theta_est               = zeros(N,1);
    theta_est(weights > 0.5)= 1;
    theta_est(weights <= 0.5)= -1;
    solution.theta_est      = theta_est;
    solution.R_est          = R_est;
    solution.t_est          = t_est;
    solution.itr            = itr;
    solution.divFactor      = divFactor;
    solution.time_gnc       = time_gnc;
    solution.f_est          = f_est;
    solution.detectedOutliers = allPoints(theta_est<0);

    % print some info
    fprintf('f_est = %g, divFactor=%g, itr=%d, time_gnc=%g[s].\n',f_est,divFactor,solution.itr,time_gnc);
    fprintf('=====================================================================================\n\n\n')
end