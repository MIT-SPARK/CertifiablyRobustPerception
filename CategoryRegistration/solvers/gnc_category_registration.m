function solution = gnc_category_registration(problem,SDP,path,varargin)
%% Solve outlier-robust category registration using GNC
%% Heng Yang, July 06, 2021

params = inputParser;
params.CaseSensitive = false;
params.addParameter('lambda',0.1, @(x) isscalar(x));
params.addParameter('refine',false, @(x) islogical(x));
params.addParameter('refineversion','v2', @(x) ischar(x));
params.parse(varargin{:});
lambda  = params.Results.lambda;
refine  = params.Results.refine;
refineversion = params.Results.refineversion;

N       = problem.N;
barc2   = 1.0;
shapes  = problem.shapes;
scene   = problem.scene;
noiseBoundSq = problem.noiseBoundSq;
cBound  = problem.cBound;
tBound  = problem.translationBound;
t0      = tic;
fprintf('\n\n===================================== GNC-TLS =======================================\n')

allPoints = 1:N;  
weights   = ones(N,1);
stopTh    = 1e-20;
maxSteps  = 1e2;
divFactor = 2;
itr       = 0;
pre_f_cost= 1e6;
cost_diff = 1e6;
breakyes  = 0;

fprintf('epsilon: %3.2e, div: %1.1f, maxiters: %3.0e.\n',stopTh,divFactor,maxSteps);
fprintf('-------------------------------------------------------------------------------------\n')
fprintf(' itr |    obj        delobj  |    mu    |   sumw   |  sumout  |   maxres   |   gap   |\n')

while (true)
    if max(abs(weights)) < 1e-6
        msg      = 'GNC encounters numerical issues, the solution is likely to be wrong.';
        breakyes = 3;
    end
    if itr == maxSteps
        breakyes = 2;
        msg      = 'Maximum iterations reached.';
    end
    if cost_diff < stopTh
        breakyes = 1;
        msg      = sprintf('GNC converged %3.2e < %3.2e.',cost_diff,stopTh);
    end
    if breakyes > 0
        fprintf('%s\n',msg);
        break
    end
    
    problem.weights     = weights;
    [R,t,c,ofinfo]      = outlier_free_category_registration(problem,path,'lambda',lambda);
    shape               = combine_shapes(shapes,c);
    residuals           = sum( (scene - R*shape - t).^2, 1) / noiseBoundSq;
    f_cost              = sum(min(residuals,barc2)) + lambda * (c'*c);
    
    maxResidual = max(residuals);
    if itr < 1
        mu = max(1 / ( 5 * maxResidual / barc2 - 1 ),1e-6); % make sure mu is positive
    end

    % update weights in closed-form
    th1 = (mu+1)/mu * barc2;
    th2 = (mu)/(mu+1) * barc2; % th1 > th2
    for i = 1:N
        if residuals(i) - th1 >= 0
            weights(i) = 0;
        elseif residuals(i) - th2 <= 0
            weights(i) = 1;
        else
            weights(i) = sqrt( barc2*mu*(mu+1)/residuals(i) ) - mu;
%             assert(weights(i)>= 0 && weights(i) <=1, 'weights calculation wrong!');
        end
    end
    cost_diff   = abs(f_cost - pre_f_cost);
    fprintf('%4d | %3.4e %3.4e | %3.2e | %3.2e | %3.2e | %3.4e | %1.1e |\n',...
            itr,f_cost,cost_diff,mu,sum(weights),N-sum(weights),maxResidual,ofinfo.gap);
    
    % increase mu and compute cost difference
    mu          = mu * divFactor;
    itr         = itr + 1;
    pre_f_cost  = f_cost;
end
R_est = R;
t_est = t;
c_est = c;

theta_est               = zeros(N,1);
theta_est(weights>0.5)  = 1;
theta_est(weights<0.5)  = -1;

% if max(max(-c_est,0)) > 0
%     c_est(c_est <= 0)       = 1e-3;
%     addpath(genpath(path.manoptpath))
%     rrPar.blk               = SDP.blk;
%     rrPar.translationBound  = problem.translationBound;
%     rrPar.cBound            = problem.cBound;
%     rrPar.N                 = problem.N;
%     rrPar.K                 = problem.K;
%    
%     [R_est,t_est]  = invert_transformation(R_est,t_est);
%     [~,~,out]      = nlp_catreg_v2(SDP.C,rrPar,R_est,t_est,c_est,theta_est);
%     [out.R,out.t]  = invert_transformation(out.R,out.t);
%     rmpath(genpath(path.manoptpath))
%     
%     R_est          = out.R;
%     t_est          = out.t;
%     c_est          = out.c;
%     theta_est      = out.theta;
%     
% end


% make sure c_est is norm bounded
nncflag = 0;
if sum(c_est <= 0) > 0
    nncflag           = 1;
    fprintf('GNC c has %d entries below 0.\n',sum(c_est <= 0));
    c_est(c_est <= 0) = 0;
end
bdcflag = 0;
if norm(c_est) >= cBound
    bdcflag           = 1;
    fprintf('GNC c has norm %g > %g.\n',norm(c_est),cBound);
    c_est = c_est / norm(c_est) * cBound;
end
% make sure t_est is norm bounded
bdtflag = 0;
if norm(t_est) >= tBound
    bdtflag     = 1;
    fprintf('GNC t has norm %g > %g.\n',norm(t_est),tBound);
    t_est = t_est / norm(t_est) * tBound;
end

% re-evaluate the cost 
shape               = combine_shapes(shapes,c_est);
residuals           = sum( (scene - R_est*shape - t_est).^2, 1) / noiseBoundSq;
f_est               = sum(min(residuals,barc2)) + lambda * (c_est'*c_est);

% if refine
%     addpath(genpath(path.manoptpath))
%     rrPar.blk               = SDP.blk;
%     rrPar.translationBound  = problem.translationBound;
%     rrPar.cBound            = problem.cBound;
%     rrPar.N                 = problem.N;
%     rrPar.K                 = problem.K;
%     switch refineversion
%     case 'v1'
%         [~,fopt,out] = nlp_catreg(SDP.C,rrPar,R_est,t_est,c_est,theta_est);
%     case 'v2'
%         [R_est,t_est] = invert_transformation(R_est,t_est);
%         [~,fopt,out]  = nlp_catreg_v2(SDP.C,rrPar,R_est,t_est,c_est,theta_est);
%         [out.R,out.t] = invert_transformation(out.R,out.t);
%     otherwise
%         error('Unknown refine version.') 
%     end
%     if fopt < f_est
%         fprintf('        MANOPT cost %3.4e < GNC cost %3.4e.\n',fopt,f_est);
%         R_est   = out.R;
%         t_est   = out.t;
%         c_est   = out.c;
%         theta_est = out.theta;
%         f_est   = fopt;
%     end
%     rmpath(genpath(path.manoptpath))
% end

time_gnc = toc(t0);

solution.type           = 'GNC-TLS';
solution.weights        = weights;
solution.theta_est      = theta_est;
solution.R_est          = R_est;
solution.t_est          = t_est;
solution.c_est          = c_est;
solution.itr            = itr;
solution.divFactor      = divFactor;
solution.time           = time_gnc;
solution.f_est          = f_est;
solution.residuals      = residuals;
solution.detectedOutliers = allPoints(theta_est<0);
solution.nncflag        = nncflag;
solution.bdcflag        = bdcflag;
solution.bdtflag        = bdtflag;

% print some info
fprintf('f_est = %g, divFactor=%g, itr=%d, time_gnc=%g[s].\n',f_est,divFactor,solution.itr,time_gnc);
fprintf('=====================================================================================\n\n\n')
end