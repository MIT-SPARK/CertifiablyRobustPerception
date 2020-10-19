function solution = tls_single_rotation_averaging_sparseSOS(problem)
    % solve TLS single rotation averaging using sparseSOS
    % Heng Yang
    % 04/14/2020
    yalmip('clear')
    % copy data
    N = problem.N;
    R_measurements = problem.R_measurements;
    noiseBoundSq = problem.noiseBoundSq;
    barc2 = 1.0;

    t0 = tic;
    % define my decision variables
    r = sdpvar(9,1);
    R = reshape(r,[3,3]);
    r1 = R(:,1); r2 = R(:,2); r3 = R(:,3);
    theta = sdpvar(N,1);
    % calculate my objective function
    residuals = {};
    for i = 1:N 
        R_i = R_measurements(:,:,i);
        residual_i = norm(R - R_i,'fro')^2;
        residuals{end+1} = residual_i / noiseBoundSq;
    end
    f_cost = 0;
    for i = 1:N 
        f_cost = f_cost + (1+theta(i))/2 * residuals{i} + (1-theta(i))/2 * barc2;
    end

    % define my feasible set
    calr = [r1'*r1 - 1 == 0;...
            r2'*r2 - 1 == 0;...
            r3'*r3 - 1 == 0;...
            r1'*r2 == 0;...
            r2'*r3 == 0;...
            r3'*r1 == 0;...
            cross(r1,r2) - r3 == 0;...
            cross(r2,r3) - r1 == 0;...
            cross(r3,r1) - r2 == 0];
    caltheta = [];
    for i = 1:N
        caltheta = [caltheta; theta(i)^2 - 1 == 0];
    end

    %%%%%%%%%%%%%%%%%%%%%%%%% Term sparsity %%%%%%%%%%%%%%%%%%%%%%%
    relaxOrder = 2;
    calr = sdpvar(calr);
    caltheta = sdpvar(caltheta);

    r_mono_2 = monolist(r,2);
    theta_mono_2 = monolist(theta,2);

    r_theta{1} = monolist([r;theta],1);
    r_theta{2} = [r_theta{1}; kron(theta,r)]; % I only generate a partial list of monomials

    M_r_theta{1} = r_theta{1} * transpose(r_theta{1});
    M_r_theta{2} = r_theta{2} * transpose(r_theta{2});

    % add constraints on moment matrices
    Fmoments = ([]);
    Fmoments = Fmoments + (M_r_theta{2}>=0); % the largest moment matrix being PSD
    PSDCon_size = size(M_r_theta{2},1);
    % localizing moments for rotation r
    for i = 1:length(calr)
        Localizer = calr(i) * theta_mono_2;
        Fmoments = Fmoments + (Localizer == 0);
    end
    % localizing moments for binary theta
    for i = 1:length(caltheta)
        Localizer = caltheta(i) * r_mono_2;
        Fmoments = Fmoments + (Localizer == 0);
    end

    % replace high-order monomials by new variables
    vars = unique( getvariables(Fmoments) );
    [~,variabletype] = yalmip('monomtable');
    nonlinears = vars(find(variabletype(vars)));
    newLinear = sdpvar(length(nonlinears),1);
    
    linears = getvariables(newLinear);
    Fnew = variablereplace(Fmoments,nonlinears,linears);
    if isa(f_cost,'sdpvar')
        f_cost_new = variablereplace(f_cost,nonlinears,linears);
    end
    
	for i = 1:length(M_r_theta)
	    if isa(M_r_theta{i},'sdpvar')
	        M_r_theta{i} = variablereplace(M_r_theta{i},nonlinears,linears);
	    end
    end

    linears = recover(linears);
    nonlinears = recover(nonlinears);
    t_relax = toc(t0); 

    % solve the SDP
    options = sdpsettings('solver', 'mosek', 'verbose', 1, 'debug',1);
    optimize(Fnew,f_cost_new,options);
    t_sdp = toc(t0) - t_relax; % this is the time used to solve the SDP

    f_sdp = value(f_cost_new);
    M_r_theta_val={};
    for i = 1:length(M_r_theta)
        M_r_theta_val{i} = value(M_r_theta{i});
    end
    
    M1 = M_r_theta_val{1};
    M2 = M_r_theta_val{2};
    
    % extract solutions
    rank_moment(1) = rank(M1,1e-3);
    rank_moment(2) = rank(M2,1e-3);
    
    [V,~] = eig(M1);
    v_monomials = V(:,end);
    v_monomials = v_monomials/v_monomials(1);
    r_est = v_monomials(2:10);
    R_est = project2SO3( reshape(r_est,[3,3]) );
    r_est = R_est(:);
    theta_raw = v_monomials(11:end);
    theta_est = theta_raw;
    theta_est(theta_est<0) = -1;
    theta_est(theta_est>0) = 1;
    f_est = replace(f_cost,[r;theta],[r_est;theta_est]);
    absDualityGap = abs(f_est - f_sdp);
    relDualityGap = absDualityGap / abs(f_est);

    % log data into solution structure
    solution.type = 'TLS';
    solution.useSparsity = true;
    solution.f_est = f_est;
    solution.f_sdp = f_sdp;
    solution.relaxOrder = relaxOrder;
    solution.absDualityGap = absDualityGap;
    solution.relDualityGap = relDualityGap;
    solution.R_est = R_est;
    solution.theta_est = theta_est;
    solution.theta_raw = theta_raw;
    solution.Moments = M_r_theta_val;
    solution.rank_moment = rank_moment;
    solution.t_relax = t_relax;
    solution.t_sdp = t_sdp;
    solution.PSDCon_size = PSDCon_size;
    allPoints = 1:N;
    solution.detectedOutliers = allPoints(theta_est<0);
    
    % print some info
    fprintf('======================== Sparse TLS Relaxation ========================\n')
    fprintf('PSD constraint size:'); fprintf('%d ',PSDCon_size); fprintf('.\n');
    fprintf('f_sdp = %g, f_est = %g, absDualityGap = %g, relDualityGap = %g, rank_moment = ',...
        f_sdp, f_est, absDualityGap, relDualityGap);
    fprintf('%g ',rank_moment);
    fprintf('.\n')
    fprintf('Relaxation time: %g[s], SDP time: %g[s].\n',t_relax,t_sdp);
    fprintf('=======================================================================\n')

end
