function solution = tls_single_rotation_averaging_gnc(problem)
    % use GNC-TLS to solve single rotation averaging problem
    % Heng Yang
    % 04/20/2020
    % copy data
    N = problem.N;
    R_measurements = problem.R_measurements;
    noiseBoundSq = problem.noiseBoundSq;
    barc2 = 1.0;

    t0 = tic;

    allPoints = 1:N;  
    weights = ones(N,1);
    stopTh = 1e-20;
    maxSteps = 1e5;
    divFactor = 1.4;
    itr = 0;
    pre_f_cost = inf;
    f_cost = inf;
    cost_diff = inf;

    while itr < maxSteps && cost_diff > stopTh
        if max(abs(weights)) < 1e-6
            disp('Weights vanish, GNC failed.')
            break;
        end
        
        % fix weights and solve for R using closed form solution
        [R_est, residuals] = weighted_single_rotation_averaging(R_measurements,weights);
        residuals = residuals / noiseBoundSq; 
        % f_cost = sum(min(residuals,barc2)); % this is the TLS cost
        f_cost = weights' * residuals;
        % fix R and residuals, and update weights in closed-form
        if itr < 1
            maxResidual = max(residuals);
            mu = max(1 / ( 5 * maxResidual / barc2 - 1 ),1e-6); % make sure mu is positive
            cprintf('Keywords', 'maxResidual=%g, set mu=%g.\n',maxResidual, mu);
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

        % increase mu and compute cost difference
        cost_diff = abs(f_cost -pre_f_cost);
        mu = mu * divFactor;
        itr = itr + 1;
        pre_f_cost = f_cost;
    end

    f_est = sum(min(residuals,barc2)); % this is the TLS cost
    R_est = R_est;
    t_gnc = toc(t0);
    
    theta_est = zeros(N,1);
    theta_est(weights > 0.5) = 1;
    theta_est(weights < 0.5) = -1;
    
    anglediff = zeros(N,1);
    for i=1:N
        Rdiff = R_est'*R_measurements(:,:,i);
        Rdiff_axang = rotm2axang(Rdiff);
        anglediff(i) = Rdiff_axang(end);
    end
    anglediff_inliers = anglediff(theta_est>0);

    solution.type = 'GNC-TLS';
    solution.nrInliers = length(anglediff_inliers);
    solution.weights = weights;
    solution.theta_est = theta_est;
    solution.R_est = R_est;
    solution.itr = itr;
    solution.divFactor = divFactor;
    solution.t_gnc = t_gnc;
    solution.f_est = f_est;
    solution.detectedOutliers = allPoints(theta_est<0);
    solution.anglediff = anglediff;
    solution.anglediff_inliers = anglediff_inliers;

    % print some info
    fprintf('============================== GNC-TLS ================================\n')
    fprintf('f_est = %g, divFactor=%g, itr=%d, t_gnc=%g[s].\n',f_est,divFactor,solution.itr,t_gnc);
    fprintf('=======================================================================\n')
end

function [R_est, residuals] = weighted_single_rotation_averaging(R_measurements,weights)
    N = zeros(3,3);
    for i = 1:length(weights)
        N = N + R_measurements(:,:,i) * weights(i);
    end
    R_est = project2SO3(N);
    residuals = zeros(length(weights),1);
    for i = 1:length(weights)
        residuals(i) = norm(R_measurements(:,:,i) - R_est, 'fro')^2;
    end
end