function solution = msac_absolute_pose_estimation(problem,SDP)
%% solve outlier-robust absolute pose estimation using MSAC
%% Depends on Matlab Computer Vision Toolbox
%% Heng Yang, June 28, 2021

N             = problem.N;
meas2D        = problem.x;
meas3D        = problem.X;
noiseBoundSq  = problem.noiseBoundSq;
pixNoiseBound = problem.pixNoiseBound;

time0       = tic;
camera      = cameraParameters('IntrinsicMatrix',eye(3));
[R_sac,t_sac,inliermask,status] ...
            = estimateWorldCameraPose(meas2D',meas3D',camera,...
                                      'MaxReprojectionError',pixNoiseBound,...
                                      'Confidence',100-1e-6,...
                                      'MaxNumTrials',1e4);                             

t_sac                   = - R_sac * t_sac';

residuals_sac = zeros(N,1);
for i = 1:N
    bearingi    = [meas2D(:,i);1];
    bearingi    = bearingi / norm(bearingi);
    pointi      = R_sac * meas3D(:,i) + t_sac;
    residuals_sac(i) = (pointi' * (eye(3) - bearingi*bearingi') * pointi)/noiseBoundSq;
end
theta_sac                       = ones(N,1);
theta_sac(residuals_sac>1.0)    = -1;
f_sac                           = sum(min(residuals_sac,1.0));

% locally refine MSAC solution using the TLS cost function
rrPar.translationBound = problem.translationBound;
rrPar.depthBound       = problem.depthBound;
rrPar.blk              = SDP.blk;
rrPar.FOV              = problem.FOV;
[res,f_est]            = nlp_ape(SDP.C,rrPar,R_sac,t_sac,theta_sac);
time_sac               = toc(time0);

R_est         = res.R;
t_est         = res.t;
theta_est     = res.theta;

solution.R_sac      = R_sac;
solution.t_sac      = t_sac;
solution.theta_sac  = theta_sac;
solution.f_sac      = f_sac;

if f_est < f_sac
    fprintf('\nLocal search found a better solution than RANSAC. %g < %g.\n',f_est,f_sac);
    solution.R_est      = R_est;
    solution.t_est      = t_est;
    solution.theta_est  = theta_est;
    solution.f_est      = f_est;
else
    fprintf('\nLocal search fails to find a better solution than RANSAC. %g > %g.\n',f_est,f_sac);
    solution.R_est      = R_sac;
    solution.t_est      = t_sac;
    solution.theta_est  = theta_sac;
    solution.f_est      = f_sac;
end
    
solution.time       = time_sac;
solution.nrInliers  = sum(solution.theta_est>0);
solution.inliermask = inliermask;

            
% print some info
fprintf('\n============================== MSAC ================================\n')
fprintf('f_est = %g, status code: %d, num outliers: %d, time_msac=%g[s].\n',...
    solution.f_est,status,N-solution.nrInliers,time_sac);
fprintf('======================================================================\n')

end