function chordal_certificate = find_chordal_sos_certificate_point_cloud_registration(problem,solution)
    fprintf('============================== Chordal SOS ==============================\n');

    % copy problem data
    yalmip('clear')
    N = problem.N;
    cloudA = problem.cloudA;
    cloudB = problem.cloudB;
    noiseBoundSq = problem.noiseBoundSq;
    barc2 = 1; % since the residuals are normalized, barc2 is always 1
    translationBound = problem.translationBound;
    % copy solution data
    f_est = solution.f_est;

    %% solve the dual SDP
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

    % use the chordal SOS 
    relaxOrder = 2;
    r_t{1} = monolist([r;t],1);
    r_t{2} = monolist([r;t],2);
    % r_theta{1} = monolist([r;theta],1);
    t_theta{2} = monolist([t;theta],2);

    r_theta_chordal = {};
    r_t_theta_chordal = {};
    X1_chordal = {};
    X_chordal = {};
    for i = 1:N 
        r_theta_chordal{end+1} = [monolist(r,1);theta(i)];
        r_t_theta_chordal{end+1} = [r_t{1};theta(i);theta(i)*r;theta(i)*t];
        X1_chordal{end+1} = sdpvar(length(r_theta_chordal{i}));
        X_chordal{end+1} = sdpvar(length(r_t_theta_chordal{i}));
    end
    fprintf('Chordal SOS: PSD Constraint Sizes: %d (X1_chordal), %d (X_Chordal).\n',length(X1_chordal{1}),length(X_chordal{1}));

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
    RHS = 0;
    for i = 1:N 
        RHS = RHS + (r_t_theta_chordal{i}' * X_chordal{i} * r_t_theta_chordal{i}) + g_t * ( (r_theta_chordal{i})' * X1_chordal{i} * r_theta_chordal{i} );
    end
    coeffs = coefficients(LHS - RHS, [r;t;theta]);
    Constraints = [coeffs == 0];
    for i = 1:N 
        Constraints = [Constraints, X_chordal{i} >= 0, X1_chordal{i}>=0];
    end

    options = sdpsettings('solver', 'mosek', 'verbose', 0, 'debug',1,'dualize',1);
    t0 = tic;
    optimize(Constraints,-gam,options);
    time_chordal_sdp = toc(t0);

    gam_est = value(gam);
    lambda_r_est = value(lambda_r);
    lambda_theta_est = value(lambda_theta);

    X_chordal_est = {};
    X1_chordal_est = {};
    X1_est = zeros(10+N,10+N);
    Xp_est = zeros(13*N+40,13*N+40);
    for i = 1:N
        X_chordal_est{end+1} = value(X_chordal{i});
        X1_chordal_est{end+1} = value(X1_chordal{i});

        idx1 = [(1:10)';10+i];
        X1_est(idx1,idx1) = X1_est(idx1,idx1) + X1_chordal_est{i};
        idx2 = [(1:13)';13+i;40+N+blkIndices(i,9);10*N+40+blkIndices(i,3)];
        Xp_est(idx2,idx2) = Xp_est(idx2,idx2) + X_chordal_est{i};
    end
    
    absDualityGap = abs(gam_est - f_est);
    relDualityGap = absDualityGap/f_est;
    fprintf('Chordal SOS: absDualityGap = %g, relDualityGap = %g, SDP time = %g[s].\n',absDualityGap,relDualityGap,time_chordal_sdp);

    % convert the solution into vectorized form
    x0_free = [lambda_r_est(:); lambda_theta_est(:)];
    nrDualVars_free = length(x0_free);
    x0_t = svec(X1_est);
    dim_svect = length(x0_t);
    x0_p = svec(Xp_est);
    dim_svecp = length(x0_p);
    x0 = [x0_free;x0_t;x0_p];

    chordal_certificate.f_est = f_est;
    chordal_certificate.gam_est = gam_est;
    chordal_certificate.absDualityGap = absDualityGap;
    chordal_certificate.relDualityGap = relDualityGap;
    chordal_certificate.lambda_r_est = lambda_r_est;
    chordal_certificate.lambda_theta_est = lambda_theta_est;
    chordal_certificate.X_chordal_est = X_chordal_est;
    chordal_certificate.X1_chordal_est = X1_chordal_est;
    chordal_certificate.X1_est = X1_est;
    chordal_certificate.Xp_est = Xp_est;
    chordal_certificate.x0 = x0;
    chordal_certificate.nrDualVars_free = nrDualVars_free;
    chordal_certificate.dim_svect = dim_svect;
    chordal_certificate.dim_svecp = dim_svecp;
    chordal_certificate.time_chordal_sdp = time_chordal_sdp;

    fprintf('=========================================================================\n');
end
