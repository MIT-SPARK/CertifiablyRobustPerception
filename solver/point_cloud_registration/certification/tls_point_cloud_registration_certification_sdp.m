function certification = tls_point_cloud_registration_certification_sdp(problem,solution)
    % Given a solution to the point cloud registation problem
    % Certify the correctness of the solution
    % Compute a sub-optimality bound
    % Directly solve the dual problem and check complementary slackness
    % Heng Yang
    % 04/06/2020
    yalmip('clear')

    % copy problem data
    N = problem.N;
    cloudA = problem.cloudA;
    cloudB = problem.cloudB;
    noiseBoundSq = problem.noiseBoundSq;
    barc2 = 1; % since the residuals are normalized, barc2 is always 1
    translationBound = problem.translationBound;
    % copy solution data
    R_est = solution.R_est;
    r_est = R_est(:);
    t_est = solution.t_est;
    theta_est = ones(N,1);
    theta_est(solution.detectedOutliers) = -1.0;
    f_est = solution.f_est;

    %% solve the feasiblity SDP
    gam = sdpvar(1);
    % define variables and constraints
    r = sdpvar(9,1);
    t = sdpvar(3,1);
    theta = sdpvar(N,1);
    R = reshape(r,[3,3]);
    r1 = R(:,1); r2 = R(:,2); r3 = R(:,3);
    g_r = [1-r1'*r1;1-r2'*r2;1-r3'*r3;...
           r1'*r2;r2'*r3;r3'*r1;...
           cross(r1,r2)-r3;cross(r2,r3)-r1;cross(r3,r1)-r2];
    g_t = translationBound^2-t'*t;
    g_theta = [];
    for i = 1:N 
        g_theta =[g_theta; 1-theta(i)^2];
    end
    % calculate residuals for each pair of measurements
    residuals = {};
    for i = 1:N
        ai = cloudA(:,i);
        bi = cloudB(:,i);
        residual = bi'*bi + ai'*ai + t'*t - 2*bi'*R*ai - 2*bi'*t + 2*t'*R*ai;
        residuals{end+1} = residual / noiseBoundSq;
    end
    % compute the objective function
    f_cost = 0;
    for i = 1:N
        fi = (1 + theta(i))/2 * residuals{i} + (1 - theta(i))/2 * barc2;
        f_cost = f_cost + fi;
    end
    
    relaxOrder = 2;
    r_t{1} = monolist([r;t],1);
    r_t{2} = monolist([r;t],2);
    r_theta{1} = monolist([r;theta],1);
    t_theta{1} = monolist([t;theta],1);
    t_theta{2} = monolist([t;theta],2);
    r_t_theta{1} = monolist([r;t;theta],1);
    r_t_theta{2} = [r_t_theta{1}; kron(r,t); kron(theta,r); kron(theta,t)]; % I only generate a partial list of monomials

    X0 = sdpvar(length(r_t_theta{2}));
    X1 = sdpvar(length(r_theta{1}));
    lambda_r = sdpvar(length(t_theta{2}), length(g_r));
    lambda_theta = sdpvar(length(r_t{2}), length(g_theta));
    r_multipliers = [];
    theta_multipliers = [];
    for i = 1:length(g_r)
        r_multipliers = [r_multipliers; lambda_r(:,i)'*t_theta{2}];
    end
    for i = 1:length(g_theta)
        theta_multipliers = [theta_multipliers; lambda_theta(:,i)'*r_t{2}];
    end


    LHS = f_cost - gam - r_multipliers' * g_r - theta_multipliers' * g_theta;
    RHS = r_t_theta{2}' * X0 * r_t_theta{2} + (r_theta{1}' * X1 * r_theta{1}) * g_t;

    coeffs = coefficients(LHS - RHS, [r;t;theta]);
    Constraints = [coeffs == 0, X0>=0, X1>=0];
    r_t_theta_2_est = replace(r_t_theta{2},[r;t;theta],[r_est;t_est;theta_est]);
    r_theta_1_est = replace(r_theta{1},[r;theta],[r_est;theta_est]);
    % complementary slackness
    % Constraints = [Constraints, X0*r_t_theta_2_est == 0, X1*r_theta_1_est == 0]; 

    %% convert the data to Sedumi standard format
%     [model,recoverymodel,diagnostics] = export(Constraints,-gam,sdpsettings('solver','sedumi'));


    options = sdpsettings('solver', 'mosek', 'verbose', 1, 'debug',1,'dualize',1);
    optimize(Constraints,-gam,options);

    gam_est = value(gam);
    X0_est = value(X0);
    X1_est = value(X1);
    lambda_r_est = value(lambda_r);
    lambda_theta_est = value(lambda_theta);
    r_multipliers_est = zeros(length(r_multipliers),1);
    for i = 1:length(r_multipliers)
        r_multipliers_est(i) = replace(r_multipliers(i), [lambda_r(:,i);t;theta],[lambda_r_est(:,i);t_est;theta_est]);
    end
    theta_multipliers_est = zeros(length(theta_multipliers),1);
    for i = 1:length(theta_multipliers)
        theta_multipliers_est(i) = replace(theta_multipliers(i), [lambda_theta(:,i);r;t],[lambda_theta_est(:,i);r_est;t_est]);
    end
    slackness_X0 = norm(X0_est * r_t_theta_2_est);
    slackness_X1 = norm(X1_est * r_theta_1_est);

    
    certification.lambda_r_est = lambda_r_est;
    certification.lambda_theta_est = lambda_theta_est;
    certification.X0_est = X0_est;
    certification.X1_est = X1_est;
    certification.slackness_X0 = slackness_X0;
    certification.slackness_X1 = slackness_X1;
    certification.gam_est = gam_est;
    certification.f_est = f_est;
    certification.suboptimality = abs(gam_est - f_est)/abs(f_est);



end