function solution = tls_mesh_registration_denseSOS(problem)
    % use denseSOS to solve TLS mesh registration
    % Heng Yang
    % 03/26/2020
    yalmip('clear')
    % copy data
    N = problem.N;
    normalM = problem.normalM;
    normalP = problem.normalP;
    pointM = problem.pointM;
    pointP = problem.pointP;
    pointNoiseBoundSq = problem.pointNoiseBoundSq;
    normalNoiseBoundSq = problem.normalNoiseBoundSq;
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
        pointM_i = pointM(:,i);
        pointP_i = pointP(:,i);
        normalM_i = normalM(:,i);
        normalP_i = normalP(:,i);
        % note the objective function is different
        residual_point_i = ( normalM_i' * (R*pointP_i - pointM_i - t) )^2;
        residual_normal_i = normalP_i'*normalP_i + normalM_i'*normalM_i - 2*normalP_i'*R'*normalM_i;
        
        residuals{end+1} = (residual_point_i/pointNoiseBoundSq + residual_normal_i/normalNoiseBoundSq)/2;
    end
    f_cost = 0;
    for i = 1:N
        fi = (1 + theta(i))/2 * residuals{i} + (1 - theta(i))/2 * barc2;
        f_cost = f_cost + fi;
    end
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
    
    calX = [calr;calt;caltheta];
    
    % derive the moment relaxation
    relaxOrder = 2;
    [Fnew,f_cost_new,momentMat] = momentmodel(calX,f_cost,relaxOrder);
    t_relax = toc(t0); % this is the time used to generate the relaxation
    % solve the resulting SDP
    options = sdpsettings('solver', 'mosek', 'verbose', 1, 'debug',1);
    optimize(Fnew,f_cost_new,options);
    t_sdp = toc(t0) - t_relax; % this is the time used to solve the SDP
    
    % extract solutions and do certification
    f_sdp = value(f_cost_new);
    momentMat_val={};
    for i = 1:length(momentMat)
        momentMat_val{i} = value(momentMat{i});
    end
    
    M1 = momentMat_val{2};
    M2 = momentMat_val{3};
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
    
    % note now I need to convert R and t back
    R_est = R_est';
    t_est = R_est * t_est;
    
    % log data into solution structure
    solution.type = 'TLS';
    solution.useSparsity = false;
    solution.f_est = f_est;
    solution.f_sdp = f_sdp;
    solution.relaxOrder = relaxOrder;
    solution.absDualityGap = absDualityGap;
    solution.relDualityGap = relDualityGap;
    solution.R_est = R_est;
    solution.t_est = t_est;
    solution.theta_est = theta_est;
    solution.theta_raw = theta_raw;
    solution.Moments = momentMat_val;
    solution.rank_moment = rank_moment;
    solution.t_relax = t_relax;
    solution.t_sdp = t_sdp;
    PSDCon_size = size(M2,1);
    solution.PSDCon_size = PSDCon_size;
    
    % print some info
    fprintf('======================== Dense TLS Relaxation ========================\n')
    fprintf('PSD constraint size:'); fprintf('%d ',PSDCon_size); fprintf('.\n');
    fprintf('f_sdp = %g, f_est = %g, absDualityGap = %g, relDualityGap = %g, rank_moment = ',...
        f_sdp, f_est, absDualityGap, relDualityGap);
    fprintf('%g ',rank_moment);
    fprintf('.\n')
    fprintf('Relaxation time: %g[s], SDP time: %g[s].\n',t_relax,t_sdp);
    fprintf('======================================================================\n')

end