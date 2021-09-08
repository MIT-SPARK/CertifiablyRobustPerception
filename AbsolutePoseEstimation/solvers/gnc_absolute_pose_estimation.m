function solution = gnc_absolute_pose_estimation(problem,varargin)
%% Solve outlier-robust absolute pose estimation using graduated non-convexity
%% Heng Yang, June 28, 2021


params = inputParser;
params.CaseSensitive = false;
params.addParameter('relaxOrder',1, @(x) isscalar(x));
params.addParameter('initPose',{}, @(x) iscell(x));
params.parse(varargin{:});

relaxOrder  = params.Results.relaxOrder;
initPose    = params.Results.initPose;

N           = problem.N;
barc2       = 1.0;    
time0       = tic;
allPoints   = 1:N;
stopTh      = 1e-20;
maxSteps    = 1e5;
divFactor   = 1.4;

if ~isfield(problem,'weights'); problem.weights = ones(N,1); end

weights     = problem.weights;
itr         = 0;
pre_f_cost  = inf;
cost_diff   = inf;
while itr < maxSteps && cost_diff > stopTh
    % fix weights and solve for R using SDP relaxation
    if max(abs(weights)) < 1e-10
        disp('Weights vanish, GNC failed.')
        break;
    end

    problem.weights     = weights;
    wls_solution        = outlier_free_absolute_pose_estimation(problem,'relaxOrder',relaxOrder);
    residuals           = wls_solution.residuals; 
    f_cost              = wls_solution.f_est;
    if itr < 1
        maxResidual = max(residuals);
        mu = max(1 / ( 5 * maxResidual / barc2 - 1 ),1e-7); % make sure mu is positive
        fprintf('maxResidual=%g, set mu=%g.\n',maxResidual, mu);
    end

    % update weights in closed-form
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
     cost_diff      = abs(f_cost -pre_f_cost);
     mu             = mu * divFactor;
     itr            = itr + 1;
     pre_f_cost     = f_cost;
end

f_est       = sum(min(residuals,barc2)); % this is the TLS cost
R_est       = wls_solution.R_est;
t_est       = wls_solution.t_est;
time_gnc    = toc(time0);

solution.type               = 'GNC-TLS';
solution.weights            = weights;
theta_est                   = zeros(N,1);
theta_est(weights > 0.5)    = 1;
theta_est(weights < 0.5)    = -1;
solution.theta_est          = theta_est;
solution.R_est              = R_est;
solution.t_est              = t_est;
solution.itr                = itr;
solution.divFactor          = divFactor;
solution.time_gnc           = time_gnc;
solution.f_est              = f_est;
solution.detectedOutliers   = allPoints(theta_est<0);

% print some info
fprintf('\n============================== GNC-TLS ================================\n')
fprintf('f_est = %g, divFactor=%g, itr=%d, time_gnc=%g[s].\n',f_est,divFactor,solution.itr,time_gnc);
fprintf('=======================================================================\n')

end