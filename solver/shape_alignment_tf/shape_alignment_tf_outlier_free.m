function solution = shape_alignment_tf_outlier_free(problem)
    % outlier-free shape alignment translation free

    yalmip('clear');
    N = problem.N;
    if ~isfield(problem,'weights')
        problem.weights = ones(N,1);
    end
    weights = problem.weights;
    z = problem.z;
    B = problem.B;
    noiseBoundSq = problem.noiseBoundSq;
    weakProjection = problem.weakProjection;
    scaleBound = problem.scaleBound;
    s_lb = scaleBound(1);
    s_ub = scaleBound(2);
    
    fprintf('SATF outlier free: scale upperbound=%g.\n',s_ub);

    r = sdpvar(6,1); 
    R = reshape(r,[2,3]); % R is actually s * weakProjection * R
    r1 = R(1,:)'; r2 = R(2,:)';
    % calculate residuals for each pair of measurements
    residuals = {};
    for i = 1:N
        zi = z(:,i);
        Bi = B(:,i);
        resVec = zi - R*Bi;
        residual = (resVec'*resVec) / noiseBoundSq;
        residuals{end+1} = weights(i) * residual;
    end
    % compute the objective function
    f_cost = 0;
    for i = 1:N
        f_cost = f_cost + residuals{i};
    end

    calr = [r1'*r1 == r2'*r2;...
            r1'*r2 == 0;...
            r'*r <= 2*s_ub^2];
    calX = calr;
    % derive the moment relaxation
    relaxOrder = 1;
    [Fnew,f_cost_new,momentMat] = momentmodel(calX,f_cost,relaxOrder);

    % solve the resulting SDP
    options = sdpsettings('solver', 'mosek', 'verbose', 0, 'debug',1);
    optimize(Fnew,f_cost_new,options);
    
    % extract solutions and do certification
    f_sdp = value(f_cost_new);
    momentMat_val={};
    for i = 1:length(momentMat)
        momentMat_val{i} = value(momentMat{i});
    end
    
    M1 = momentMat_val{2};
    rank_moment(1) = rank(M1,1e-3);
    
    [V,~] = eig(M1);
    v_monomials = V(:,end);
    v_monomials = v_monomials/v_monomials(1);
    
    r_est = v_monomials(2:7);
    [s_est,R_est] = get_s_and_R(r_est);
    r_est = s_est * reshape(R_est(1:2,:),[6,1]);

    f_est = replace(f_cost,r,r_est);
    absDualityGap = abs(f_est - f_sdp);
    relDualityGap = absDualityGap / abs(f_est);

    solution.type = 'outlier-free SOS';
    solution.s_est = s_est;
    solution.R_est = R_est;
    solution.t_est = zeros(2,1); % just for completeness
    solution.f_sdp = f_sdp;
    solution.f_est = f_est;
    solution.absDualityGap = absDualityGap;
    solution.relDualityGap = relDualityGap;
    solution.rank_moment = rank_moment;

    residuals = zeros(N,1);
    for i = 1:N 
        res = norm( z(:,i) - s_est * weakProjection * R_est * B(:,i) )^2 / noiseBoundSq;
        residuals(i) = res;
    end

    solution.residuals = residuals;
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