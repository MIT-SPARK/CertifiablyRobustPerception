function solution = tls_point_cloud_registration_sparseSOS(problem)
    % use sparse SOS relaxation to solve point cloud registration
    % with TLS 
    % Heng Yang
    % 03/23/2020
    yalmip('clear')
    % copy data
    N = problem.N;
    cloudA = problem.cloudA;
    cloudB = problem.cloudB;
    noiseBoundSq = problem.noiseBoundSq;
    barc2 = 1; % since the residuals are normalized, barc2 is always 1
    translationBound = problem.translationBound;
    
    % define decision variables
    % NOTE: I tried using unit quaternion to represent rotation, but it
    % failed, that's why I switched to rotation matrix -- HY
    t0 = tic;
    r = sdpvar(9,1);
    R = reshape(r,[3,3]);
    r1 = R(:,1); r2 = R(:,2); r3 = R(:,3);
    t = sdpvar(3,1);
    theta = sdpvar(N,1);
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
        
    calt = [t'*t <= translationBound^2];
    caltheta = [];
    for i = 1:N
        caltheta = [caltheta; theta(i)^2 - 1 == 0];
    end
    %%%%%%%%%%%%%%%%%%%%%% Term Sparsity %%%%%%%%%%%%%%%%%%%%%%%%%%
    calr = sdpvar(calr);
    calt = sdpvar(calt);
    caltheta = sdpvar(caltheta);
    
    % generate monomials and moment matrices
    relaxOrder = 2;
    % r_t{1} = monolist([r;t],1);
    r_theta{1} = monolist([r;theta],1);
    % t_theta{1} = monolist([t;theta],1);
    
    % M_r_t{1} = r_t{1} * transpose(r_t{1}); % moment matrix of only r and t
    vec_M_r_t{1} = monolist([r;t],2); % vectorization of the moment matrix 
    M_r_theta{1} = r_theta{1} * transpose(r_theta{1});
    % M_t_theta{1} = t_theta{1} * transpose(t_theta{1});
    vec_M_t_theta{1} = monolist([t;theta],2); % this is the vectorization of the moment matrix;
    
    r_t_theta{1} = monolist([r;t;theta],1);
    r_t_theta{2} = [r_t_theta{1}; kron(r,t); kron(theta,r); kron(theta,t)]; % I only generate a partial list of monomials
%     r_t_theta{2} = [r_t_theta{1}; kron(theta,r); kron(theta,t)]; % I only generate a partial list of monomials
%     r_t_theta{2} = [r_t{1}; kron(theta,r); kron(theta,t)]; % I only generate a partial list of monomials
    
    M_r_t_theta{1} = r_t_theta{1} * transpose(r_t_theta{1});
    M_r_t_theta{2} = r_t_theta{2} * transpose(r_t_theta{2});
    
    % add constraints on moment matrices
    Fmoments = ([]);
    Fmoments = Fmoments + (M_r_t_theta{2}>=0); % the largest moment matrix being PSD
    PSDCon_size = size(M_r_t_theta{2},1);
    % localizing moments for rotation r
    for i = 1:length(calr)
        Localizer = calr(i) * vec_M_t_theta{1};
        Fmoments = Fmoments + (Localizer == 0);
    end
    % localizing moments for translation t
    for i = 1:length(calt)
        Localizer = calt(i) * M_r_theta{1};
        Fmoments = Fmoments + (Localizer >= 0);
        PSDCon_size = [PSDCon_size, size(Localizer,1)];
    end
    % localizing moments for binary theta
    for i = 1:length(caltheta)
        Localizer = caltheta(i) * vec_M_r_t{1};
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
    
	for i = 1:length(M_r_t_theta)
	    if isa(M_r_t_theta{i},'sdpvar')
	        M_r_t_theta{i} = variablereplace(M_r_t_theta{i},nonlinears,linears);
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
    M_r_t_theta_val={};
    for i = 1:length(M_r_t_theta)
        M_r_t_theta_val{i} = value(M_r_t_theta{i});
    end
    
    M1 = M_r_t_theta_val{1};
    M2 = M_r_t_theta_val{2};
    
    % extract solutions
    rank_moment(1) = rank(M1,1e-3);
    rank_moment(2) = rank(M2,1e-3);
    
    [V,~] = eig(M1);
    v_monomials = V(:,end);
    v_monomials = v_monomials/v_monomials(1);
    r_est = v_monomials(2:10);
    R_est = project2SO3( reshape(r_est,[3,3]) );
    r_est = R_est(:);
    t_est = v_monomials(11:13);
    theta_raw = v_monomials(14:end);
    theta_est = theta_raw;
    theta_est(theta_est<0) = -1;
    theta_est(theta_est>0) = 1;
    f_est = replace(f_cost,[r;t;theta],[r_est;t_est;theta_est]);
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
    solution.t_est = t_est;
    solution.theta_est = theta_est;
    solution.theta_raw = theta_raw;
    solution.Moments = M_r_t_theta_val;
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

    
    
    