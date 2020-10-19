function solution = tls_shape_alignment_tf_denseSOS(problem)
    % use denseSOS to solve TLS shape alignment translation free
    % Heng Yang
    % 05/16/2020
    yalmip('clear')
    
    % copy data
    N = problem.N;
    z = problem.z;
    B = problem.B;
    noiseBoundSq = problem.noiseBoundSq;
    weakProjection = problem.weakProjection;
    barc2 = 1;
    
    
    scaleBound = problem.scaleBound;
    s_ub = scaleBound(2);
    fprintf('Upper bound of scale in SOS relaxation translation free: %g.\n',s_ub);
    
    % define decision variables
    t0 = tic;
    r = sdpvar(6,1); 
    R = reshape(r,[2,3]);
    r1 = R(1,:)'; r2 = R(2,:)';
    theta = sdpvar(N,1);
    % calculate residuals for each pair of measurements
    residuals = {};
    for i = 1:N
        zi = z(:,i);
        Bi = B(:,i);
        resVec = zi - R*Bi;
        residual = resVec' * resVec;
        residuals{end+1} = residual / noiseBoundSq;
    end
    % compute the objective function
    f_cost = 0;
    for i = 1:N
        fi = (1 + theta(i))/2 * residuals{i} + (1 - theta(i))/2 * barc2;
        f_cost = f_cost + fi;
    end
    % define my feasible set
    calr = [r1'*r1 == r2'*r2;...
            r1'*r2 == 0;...
            r'*r <= 2*s_ub^2];
    caltheta = [];
    for i = 1:N
        caltheta = [caltheta; theta(i)^2 == 1];
    end

    calX = [calr;caltheta];
    


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
    
    r_est = v_monomials(2:7);
    [s_est,R_est] = get_s_and_R(r_est);
    r_est = s_est * reshape(R_est(1:2,:),[6,1]);
    
    theta_raw = v_monomials(8:end);
    theta_est = theta_raw;
    theta_est(theta_est<0) = -1;
    theta_est(theta_est>0) = 1;
    
    f_est = replace(f_cost,[r;theta],[r_est;theta_est]);
    absDualityGap = abs(f_est - f_sdp);
    relDualityGap = absDualityGap / abs(f_est);
    
    solution.type = 'TLS';
    solution.useSparsity = false;
    solution.f_sdp = f_sdp;
    solution.f_est = f_est;
    solution.relaxOrder = relaxOrder;
    solution.absDualityGap = absDualityGap;
    solution.relDualityGap = relDualityGap;
    solution.s_est = s_est;
    solution.R_est = R_est;
    solution.t_est = zeros(2,1);
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


function [s,R] = get_s_and_R(r)
    R = reshape(r,[2,3]);
    [U,S,V] = svd(R);
    sigma_1 = S(1,1); sigma_2 = S(2,2);
    s = (sigma_1 + sigma_2) / 2;
    S_bar = [s, 0, 0; 0, s, 0];
    R = U * S_bar * V';
    R = R / s;
    % populate the last row
    R = [R; cross(R(1,:),R(2,:))];
end
