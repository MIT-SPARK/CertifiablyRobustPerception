function solution = tls_shape_alignment_tf_sparseSOS(problem)
    % use sparseSOS to solve TLS shape alignment
    % Heng Yang
    % 04/26/2020
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
    
    %%%%%%%%%%%%%%%%%%%%% derive sparse relaxation %%%%%%%%%%%%%%%%%%%%%%%%
    calr = sdpvar(calr);
    caltheta = sdpvar(caltheta);
    
    iseq_calr = [ones(2,1);0];
    
    % generate monomials and moment matrices
    relaxOrder = 2;
    % multipliers for calr
    theta_mono{2} = monolist(theta,2);
    theta_mono{1} = monolist(theta,1);
    
    M_theta{1} = theta_mono{1} * theta_mono{1}';
    
    % multipliers for caltheta
    r_mono{2} = monolist(r,2);
    
    % multipliers for 1
    r_theta{1} = monolist([r;theta],1);
    
    % even more basis reduction
%     r_theta{2} = [r_theta{1};kron(theta,r)];
    
    % standard basis reduction as written in paper
    r_mono_2 = r_mono{2}(1+6+1:end);
    r_theta{2} = [r_theta{1};r_mono_2;kron(theta,r)];
    
    M_r_theta{1} = r_theta{1} * r_theta{1}';
    M_r_theta{2} = r_theta{2} * r_theta{2}';
    
    r_theta_full_mono_2 = monolist([r;theta],2);
    
    
    % add constraints on moment matrices
    Fmoments = ([]);
    Fmoments = Fmoments + (M_r_theta{2}>=0); % the largest moment matrix being PSD
    PSDCon_size = size(M_r_theta{2},1);
    % localizing moments for rotation r
    for i = 1:length(calr)
        if iseq_calr(i)
            % Even more basis reduction
%             Localizer = calr(i) * theta_mono{2};
            
            % standard basis reduction as written in paper
            Localizer = calr(i) * r_theta_full_mono_2;
            
            Fmoments = Fmoments + (Localizer == 0);
        else
            % Even more basis reduction
%             Localizer = calr(i) * M_theta{1};
            
            % standard basis reduction as written in paper
            Localizer = calr(i) * M_r_theta{1};
            
            Fmoments = Fmoments + (Localizer >= 0);
            PSDCon_size = [PSDCon_size, size(Localizer,1)];
        end        
    end
    
    % localizing moments for binary theta
    for i = 1:length(caltheta)
        Localizer = caltheta(i) * r_mono{2};
        
%         Localizer = caltheta(i) * r_theta_full_mono_2;
        
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
    %%%%%%%%%%%%%%%%%%%%%%%% end derivation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    solution.useSparsity = true;
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
    solution.Moments = M_r_theta_val;
    solution.rank_moment = rank_moment;
    solution.t_relax = t_relax;
    solution.t_sdp = t_sdp;
    solution.PSDCon_size = PSDCon_size;
    
    % print some info
    fprintf('======================== Sparse TLS Relaxation ========================\n')
    fprintf('PSD constraint size:'); fprintf('%d ',PSDCon_size); fprintf('.\n');
    fprintf('f_sdp = %g, f_est = %g, absDualityGap = %g, relDualityGap = %g, rank_moment = ',...
        f_sdp, f_est, absDualityGap, relDualityGap);
    fprintf('%g ',rank_moment);
    fprintf('.\n')
    fprintf('Relaxation time: %g[s], SDP time: %g[s].\n',t_relax,t_sdp);
    fprintf('======================================================================\n')
end

% function R = get_R(r)
%     R2 = reshape(r,[2,3]);
%     [U,S,V] = svd(R2);
%     S = [1,0,0;0,1,0];
%     R2 = U * S * V';
%     % populate the last row
%     R = [R2; cross(R2(1,:),R2(2,:))];
% end

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