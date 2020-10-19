function chordal_certificate = find_chordal_sos_certificate_single_rotation_averaging(problem,solution)

    fprintf('============================== Chordal SOS ==============================\n');
    % copy problem data
    yalmip('clear')
    N = problem.N;
    R_measurements = problem.R_measurements;
    noiseBoundSq = problem.noiseBoundSq;
    barc2 = 1; % since the residuals are normalized, barc2 is always 1
    % copy solution data
    f_est = solution.f_est;
    R_est = solution.R_est;
    r_est = R_est(:);
    theta_est = solution.theta_est;
    vp_est = [1;r_est;theta_est;kron(theta_est,r_est)];

    %% solve the chordal sparse version of the dual SDP
    gam = sdpvar(1);
    % define variables and constraints
    r = sdpvar(9,1);
    theta = sdpvar(N,1);
    R = reshape(r,[3,3]);
    r1 = R(:,1); r2 = R(:,2); r3 = R(:,3);
    g_r = [1-r1'*r1;1-r2'*r2;1-r3'*r3;...
           r1'*r2;r2'*r3;r3'*r1;...
           cross(r1,r2)-r3;cross(r2,r3)-r1;cross(r3,r1)-r2];
    g_theta = [];
    for i = 1:N 
        g_theta =[g_theta; 1-theta(i)^2];
    end

    % calculate the residuals polynomials
    residuals = {};
    for i = 1:N 
        ri = reshape(R_measurements(:,:,i),[9,1]);
        residuals{end+1} = (r - ri)' * (r-ri) / noiseBoundSq;
    end
    f_cost = 0;
    for i = 1:N 
        f_cost = f_cost + (1+theta(i))/2 * residuals{i} + (1-theta(i))/2 * barc2;
    end

    % use the chordal SOS 
    relaxOrder = 2;
    r_mono{1} = monolist(r,1);
    r_mono{2} = monolist(r,2);
    theta_mono_2 = monolist(theta,2);
    r_theta_chordal = {};
    X_chordal = {};
    for i = 1:N 
        r_theta_chordal{end+1} = [r_mono{1};theta(i);theta(i)*r];
        X_chordal{end+1} = sdpvar(length(r_theta_chordal{i}));
    end
    fprintf('Chordal SOS: PSD Constraint Size: %d (X_Chordal).\n',length(X_chordal{1}));
    lambda_r = sdpvar(length(theta_mono_2), length(g_r));
    lambda_theta = sdpvar(length(r_mono{2}), length(g_theta));
    r_multipliers = [];
    theta_multipliers = [];
    for i = 1:length(g_r)
        r_multipliers = [r_multipliers; lambda_r(:,i)'*theta_mono_2];
    end
    for i = 1:length(g_theta)
        theta_multipliers = [theta_multipliers; lambda_theta(:,i)'*r_mono{2}];
    end
    LHS = f_cost - gam - r_multipliers' * g_r - theta_multipliers' * g_theta;
    RHS = 0;
    for i = 1:N 
        RHS = RHS + (r_theta_chordal{i}' * X_chordal{i} * r_theta_chordal{i});
    end
    coeffs = coefficients(LHS - RHS, [r;theta]);
    Constraints = [coeffs == 0];
    for i = 1:N 
        Constraints = [Constraints, X_chordal{i} >= 0];
    end

    options = sdpsettings('solver', 'mosek', 'verbose', 0, 'debug',1,'dualize',1);
    t0 = tic;
    optimize(Constraints,-gam,options);
    time_chordal_sdp = toc(t0);

    gam_est = value(gam);
    lambda_r_est = value(lambda_r);
    lambda_theta_est = value(lambda_theta);

    X_chordal_est = {};
    Xp_est = zeros(10*N+10,10*N+10);
    for i = 1:N
        X_chordal_est{end+1} = value(X_chordal{i});

        idx2 = [(1:10)';10+i;10+N+blkIndices(i,9)];
        Xp_est(idx2,idx2) = Xp_est(idx2,idx2) + X_chordal_est{i};
    end

    absDualityGap = abs(gam_est - f_est);
    relDualityGap = absDualityGap/f_est;
    PDSlackness = norm( Xp_est * vp_est );
    fprintf('Chordal SOS: absDualityGap = %g, relDualityGap = %g, PDSlackness = %g, SDP time = %g[s].\n',absDualityGap,relDualityGap,PDSlackness,time_chordal_sdp);

    % convert the solution into vectorized form
    x0_free = [lambda_r_est(:); lambda_theta_est(:)];
    nrDualVars_free = length(x0_free);
    x0_p = svec(Xp_est);
    dim_svecp = length(x0_p);
    x0 = [x0_free;x0_p];

    chordal_certificate.vp_est = vp_est;
    chordal_certificate.f_est = f_est;
    chordal_certificate.gam_est = gam_est;
    chordal_certificate.absDualityGap = absDualityGap;
    chordal_certificate.relDualityGap = relDualityGap;
    chordal_certificate.PDSlackness = PDSlackness;
    chordal_certificate.lambda_r_est = lambda_r_est;
    chordal_certificate.lambda_theta_est = lambda_theta_est;
    chordal_certificate.X_chordal_est = X_chordal_est;
    chordal_certificate.Xp_est = Xp_est;
    chordal_certificate.x0 = x0;
    chordal_certificate.nrDualVars_free = nrDualVars_free;
    chordal_certificate.dim_svecp = dim_svecp;
    chordal_certificate.time_chordal_sdp = time_chordal_sdp;

    fprintf('=========================================================================\n');
end