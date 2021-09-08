function solution = gnc_single_rotation_averaging(problem)
    % use GNC-TLS to solve single rotation averaging problem
    % Heng Yang
    % 07/17/2021
    
    % copy data
    N               = problem.N;
    R_measurements  = problem.R_measurements;
    noiseBoundSq    = problem.noiseBoundSq;
    barc2           = 1.0;

    t0 = tic;
    fprintf('\n\n===================================== GNC-TLS =======================================\n')

    allPoints   = 1:N;  
    weights     = ones(N,1);
    stopTh      = 1e-20;
    maxSteps    = 1e5;
    divFactor   = 1.4;
    itr         = 0;
    pre_f_cost  = 1e6;
    f_cost      = 1e6;
    cost_diff   = 1e6;
    R_est       = eye(3);
    
    fprintf('epsilon: %3.2e, div: %1.1f, maxiters: %3.0e.\n',stopTh,divFactor,maxSteps);
    fprintf('-------------------------------------------------------------------------------------\n')
    fprintf(' itr |    obj        delobj  |    mu    |   sumw   |  sumout  |   maxres   |   gap   |\n')
    
    breakyes    = 0;
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
        
        % fix weights and solve for R using closed form solution
        [R_est, residuals] = weighted_single_rotation_averaging(R_measurements,weights);
        residuals          = residuals / noiseBoundSq; 
        % f_cost = sum(min(residuals,barc2)); % this is the TLS cost
        f_cost             = weights' * residuals;
        f_cost_TLS         = sum(min(residuals,barc2));
        % fix R and residuals, and update weights in closed-form
        maxResidual        = max(residuals);
        if itr < 1
            mu = max(1 / ( 5 * maxResidual / barc2 - 1 ),1e-6); % make sure mu is positive
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
        cost_diff  = abs(f_cost - pre_f_cost);
        
        fprintf('%4d | %3.4e %3.4e | %3.2e | %3.2e | %3.2e | %3.4e | %1.1e |\n',...
            itr,f_cost_TLS,cost_diff,mu,sum(weights),N-sum(weights),maxResidual,0);
        
        % increase mu and compute cost difference
        mu         = mu * divFactor;
        itr        = itr + 1;
        pre_f_cost = f_cost;
    end

    f_est    = sum(min(residuals,barc2)); % this is the TLS cost
    time_gnc = toc(t0);
    
    theta_est = zeros(N,1);
    theta_est(weights > 0.5) = 1;
    theta_est(weights < 0.5) = -1;

    solution.type = 'GNC-TLS';
    solution.nrInliers = sum(weights > 0.5);
    solution.weights = weights;
    solution.theta_est = theta_est;
    solution.R_est = R_est;
    solution.itr = itr;
    solution.divFactor = divFactor;
    solution.time_gnc = time_gnc;
    solution.f_est = f_est;
    solution.detectedOutliers = allPoints(theta_est<0);

    % print some info
    fprintf('f_est = %g, divFactor=%g, itr=%d, time_gnc=%g[s].\n',f_est,divFactor,solution.itr,time_gnc);
    fprintf('=====================================================================================\n\n\n')
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