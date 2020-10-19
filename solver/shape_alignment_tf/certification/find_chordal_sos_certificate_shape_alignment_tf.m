function chordal_certificate = find_chordal_sos_certificate_shape_alignment_tf(problem,solution)
    % shape alignment translation free chordal certificate
    fprintf('============================== Chordal SOS ==============================\n');
    yalmip('clear');
    % copy problem data
    N = problem.N;
    z = problem.z;
    B = problem.B;
    noiseBoundSq = problem.noiseBoundSq;
    s_ub = problem.scaleBound(2);

    barc2 = 1.0;
    % copy solution data
    f_est = solution.f_est;

     %% solve the chordal sparse version of the dual SDP
     gam = sdpvar(1);
     % define variables and constraints
     r = sdpvar(6,1);
     theta = sdpvar(N,1);
     R = reshape(r,[2,3]);
     r1 = R(1,:)'; r2 = R(2,:)';
     g_r_eq = [r1'*r1-r2'*r2;...
               r1'*r2];

    g_r_ineq = [2*s_ub^2 - r'*r];

     g_theta = [];
     for i = 1:N 
         g_theta =[g_theta; 1-theta(i)^2];
     end

    % calculate the residuals polynomials
    residuals = {};
    for i = 1:N 
        zi = z(:,i);
        Bi = B(:,i);
        resVec = zi - R*Bi;
        residuals{end+1} = (resVec' * resVec) / noiseBoundSq;
    end
    f_cost = 0;
    for i = 1:N 
        f_cost = f_cost + (1+theta(i))/2 * residuals{i} + (1-theta(i))/2 * barc2;
    end

    % use the chordal SOS 
    relaxOrder = 2;
    r_mono{1} = monolist(r,1);
    r_mono{2} = monolist(r,2);
    theta_mono{2} = monolist(theta,2);

    theta_chordal = {};
    r_theta_chordal = {};
    Xr_chordal = {};
    Xp_chordal = {};
    for i = 1:N
        theta_chordal{end+1} = [1;theta(i)]; % (1; theta_i) * (1; theta_i)'
        r_theta_chordal{end+1} = [1;r;theta(i);theta(i)*r]; % (1; r; theta_i; theta_i*r) * (1; r; theta_i; theta_i*r)'

        Xr_chordal{end+1} = sdpvar(length(theta_chordal{i}));
        Xp_chordal{end+1} = sdpvar(length(r_theta_chordal{i}));
    end
    fprintf('Chordal SOS: PSD Constraint Size: %d (Xr_chordal), %d (Xp_Chordal).\n',length(Xr_chordal{1}),length(Xp_chordal{1}));

    lambda_r = sdpvar(length(theta_mono{2}), length(g_r_eq));
    lambda_theta = sdpvar(length(r_mono{2}), length(g_theta));
    r_multipliers = [];
    theta_multipliers = [];
    for i = 1:length(g_r_eq)
        r_multipliers = [r_multipliers; lambda_r(:,i)'*theta_mono{2}];
    end
    for i = 1:length(g_theta)
        theta_multipliers = [theta_multipliers; lambda_theta(:,i)'*r_mono{2}];
    end

    LHS = f_cost - gam - r_multipliers' * g_r_eq - theta_multipliers' * g_theta;
    RHS = 0;
    for i = 1:N 
        RHS = RHS + (r_theta_chordal{i}' * Xp_chordal{i} * r_theta_chordal{i}) + ...
                     g_r_ineq * (theta_chordal{i}' * Xr_chordal{i} * theta_chordal{i});
    end
    coeffs = coefficients(LHS - RHS, [r;theta]);
    Constraints = [coeffs == 0];
    for i = 1:N 
        Constraints = [Constraints, Xp_chordal{i} >= 0, Xr_chordal{i} >= 0];
    end

    options = sdpsettings('solver', 'mosek', 'verbose', 0, 'debug',1,'dualize',1);
    t0 = tic;
    optimize(Constraints,-gam,options);
    time_chordal_sdp = toc(t0);

    gam_est = value(gam);
    lambda_r_est = value(lambda_r);
    lambda_theta_est = value(lambda_theta);

    Xp_chordal_est = {};
    Xp_est = zeros(7*N+7,7*N+7); % 1; r; theta; kron(theta, r), 1+6+N+6N = 7N + 7
    Xr_chordal_est = {}; 
    Xr_est = zeros(1+N,1+N); % 1; theta, 1+N
    for i = 1:N
        Xp_chordal_est{end+1} = value(Xp_chordal{i});
        Xr_chordal_est{end+1} = value(Xr_chordal{i});

        idx1 = [1;1+i];
        Xr_est(idx1,idx1) = Xr_est(idx1,idx1) + Xr_chordal_est{i};

        idx2 = [(1:7)';7+i;7+N+blkIndices(i,6)];
        Xp_est(idx2,idx2) = Xp_est(idx2,idx2) + Xp_chordal_est{i};
    end

    absDualityGap = abs(gam_est - f_est);
    relDualityGap = absDualityGap/f_est;
    fprintf('Chordal SOS: absDualityGap = %g, relDualityGap = %g, SDP time = %g[s].\n',absDualityGap,relDualityGap,time_chordal_sdp);

    % convert the solution into vectorized form
    x0_free = [lambda_r_est(:); lambda_theta_est(:)];
    nrDualVars_free = length(x0_free);
    x0_r = svec(Xr_est);
    dim_svecr = length(x0_r);
    x0_p = svec(Xp_est);
    dim_svecp = length(x0_p);
    x0 = [x0_free;x0_r;x0_p];

    chordal_certificate.f_est = f_est;
    chordal_certificate.gam_est = gam_est;
    chordal_certificate.absDualityGap = absDualityGap;
    chordal_certificate.relDualityGap = relDualityGap;
    chordal_certificate.lambda_r_est = lambda_r_est;
    chordal_certificate.lambda_theta_est = lambda_theta_est;
    chordal_certificate.Xp_est = Xp_est;
    chordal_certificate.Xr_est = Xr_est;
    chordal_certificate.x0 = x0;
    chordal_certificate.nrDualVars_free = nrDualVars_free;
    chordal_certificate.dim_svecp = dim_svecp;
    chordal_certificate.dim_svecr = dim_svecr;
    chordal_certificate.time_chordal_sdp = time_chordal_sdp;

    fprintf('=========================================================================\n');







end