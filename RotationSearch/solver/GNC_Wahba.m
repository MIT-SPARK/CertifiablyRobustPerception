function [R_gnc, info_gnc] = GNC_Wahba(v1, v2, barc2, div_factor, varargin)
% solve the robust Wahba problem using graduated non-convexity

if size(v1,1) ~= 3 || size(v2,1) ~= 3
    error('GNC_Wahba works on N 3D vectors (3xN matrices)')
end

params = inputParser;
params.CaseSensitive = false;
params.addParameter('plotCostTraj',false, @(x) islogical(x));
params.parse(varargin{:});

plotCostTraj = params.Results.plotCostTraj;

t_start = tic;

nrPoints = size(v1,2);
allPoints = 1:nrPoints;

weights = ones(1,nrPoints);
stopTh = 1e-16;
maxSteps = 1e2;
itr = 0;
pre_TLS_cost = inf;
cost_diff = inf;

TLS_cost_traj = [];

while itr < maxSteps && cost_diff > stopTh
    % fix weights and solve for R using procrustes
    weighted_model = zeros(size(v1));
    weighted_scene = weighted_model;
    for i = 1:nrPoints
        weighted_model(:,i) = weights(i)^(0.5) * v1(:,i);
        weighted_scene(:,i) = weights(i)^(0.5) * v2(:,i);
    end
    % do SVD of data matrix
    R_est_cf = Wahba_closed_form(weighted_model, weighted_scene);
    
    % fix R, solve for weights in closed form
    residuals = zeros(1,nrPoints);
    for i = 1:nrPoints
        residuals(i) = norm( v2(:,i) - R_est_cf * v1(:,i) )^2;
    end

    if itr < 1
        maxResidual = max(residuals);
        mu = max(1 / ( 5 * maxResidual / barc2 - 1 ), 1e-3 );
%         fprintf('GNC-TLS first iteration: maxResidual=%g, set mu=%g.\n',maxResidual, mu);
    end

    th1 = (mu+1)/mu * barc2;
    th2 = (mu)/(mu+1) * barc2; % th1 > th2
    for i = 1:nrPoints
        if residuals(i) - th1 >= 0
            weights(i) = 0;
        elseif residuals(i) - th2 <= 0
            weights(i) = 1;
        else
            weights(i) = sqrt( barc2*mu*(mu+1)/residuals(i) ) - mu;
            assert(weights(i)>= 0 && weights(i) <=1, 'weights calculation wrong!');
        end
    end

    TLS_cost = weights * residuals' + ( nrPoints - sum(weights) ) * barc2;

    TLS_cost_traj = [TLS_cost_traj, TLS_cost];

    cost_diff = abs(TLS_cost - pre_TLS_cost);

    % increase mu
    mu = mu * div_factor;
    itr = itr + 1;
    pre_TLS_cost = TLS_cost;
end

R_gnc = R_est_cf;

detected_outliers = allPoints(weights < 0.5);
detected_inliers = allPoints(weights > 0.5);

theta_gnc = ones(nrPoints,1);
theta_gnc(detected_outliers) = -1;

t_gnc = toc(t_start);

info_gnc.maxSteps = maxSteps;
info_gnc.weights = weights;
info_gnc.itr = itr;
info_gnc.finalMu = mu;
info_gnc.TLS_cost_traj = TLS_cost_traj;
info_gnc.detected_outliers = detected_outliers;
info_gnc.detected_inliers = detected_inliers;
info_gnc.weights_sum = sum(weights);
info_gnc.theta_gnc = theta_gnc;
info_gnc.t_gnc = t_gnc;

if plotCostTraj
    figure
    set(gcf,'position',[400,400,400,600])
    subplot(2,1,1)
    plot(weights,'linewidth',2);
    ylim([0,1])
    title('weights returned by GNC-TLS')
    xlabel('Index')
    ylabel('Weight')
    
    subplot(2,1,2)
    plot(TLS_cost_traj,'linewidth',2);
    title('Trajectory of TLS cost')
    xlabel('Iteration')
    ylabel('Cost')
end
