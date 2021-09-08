function problem = gen_absolute_pose_estimation(problem)
%% generate a random absolute pose estimation problem
%% Also known as the perspective-n-points (PnP) problem
%% 3D points are centered at zero (reasonable for a CAD model)
%% Translation bounded, and depth in front of camera
%% Heng Yang
%% June 26, 2021

if ~isfield(problem, 'N'); error('Please use problem.N to specify the number of correspondences.'); end
if ~isfield(problem, 'FOV'); problem.FOV = 90; end % default field of view 90 degrees
if ~isfield(problem, 'noiseSigma'); problem.noiseSigma = 0.01; end
if ~isfield(problem, 'translationBound'); problem.translationBound = 20.0; end % maximum translation 20 meters
if ~isfield(problem, 'outlierRatio'); problem.outlierRatio = 0.0; end
if ~isfield(problem, 'minDepth'); problem.minDepth = 1.0; end % minimum depth of the 3D points, not too close to the camera

N                 = problem.N;
FOV               = problem.FOV;
noiseSigma        = problem.noiseSigma;
translationBound  = problem.translationBound;
outlierRatio      = problem.outlierRatio;
minDepth          = problem.minDepth;

% first generate a random 3D model, centered at zero
X   = randn(3,N);
X   = X - mean(X,2); % centered at zero

% generate random rotation
R       = rand_rotation;
% generate a random translation direction in the FOV cone
halfFOV = deg2rad(FOV/2);
tmp     = tan(rand*halfFOV);
alpha   = 2*pi*rand;
t       = [tmp*cos(alpha);tmp*sin(alpha);1];
t       = t/norm(t); % unit vector as translation direction

% apply R and t such that X_c is in front of the camera
RX      = R*X;
RX3     = min(RX(3,:));
tscale  = translationBound * rand;
config  = true;
while config
    if tscale * t(end) + RX3 > minDepth
        config = false;
    else
        tscale = tscale * 1.1; % increase tscale if not large enough
        if tscale > translationBound
            error('Translation bound too small.');
        end
    end
end
t       = tscale * t; % scaled real translation
X_c     = R*X + t; % 3D points in camera frame

assert(t(1)^2 + t(2)^2 <= t(3)^2, 'Center of object not in FOV cone.');

% generate projections
x_h     = X_c ./ X_c(3,:);
x       = x_h(1:2,:); % 2 by N

% add noise
x       = x + noiseSigma * randn(2,N);

% add outliers
nrOutliers = round(N*outlierRatio);
if nrOutliers > 0
    fprintf('Absolute Pose Estimation: generate %d random outliers.\n',nrOutliers);
    % random outliers in the FOV cone
    tmp             = tan( rand(nrOutliers,1) * halfFOV );
    alpha           = 2*pi*rand(nrOutliers,1);
    outliers        = [tmp.*cos(alpha),tmp.*sin(alpha)];
    outlierIDs      = (N-nrOutliers+1):N;
    x(:,outlierIDs) = outliers';
else
    outlierIDs = [];
end

% compute noiseBound, noiseBound should not be too small, otherwise solvers
% may encounter difficulties
pixNoiseBound     = sqrt(noiseSigma^2 * chi2inv(0.9999,2));
noiseBoundSq      = max(0.2^2,(noiseSigma^2 * chi2inv(0.99,2)));
noiseBound        = sqrt(noiseBoundSq);
fprintf('N: %d, outlierRatio: %g, noiseBound: %g, noiseBoundSq: %g, pixNoiseBound: %g, translationBound: %g.\n',...
    N,outlierRatio,noiseBound,noiseBoundSq,pixNoiseBound,translationBound);

% compute residuals
residuals       = compute_point_line_distances(X_c,x);

R_gt            = R;
t_gt            = t;
problem.X       = X; % 3D points
problem.x       = x; % noisy and outlier-corrupted 2D measurements
problem.R_gt    = R_gt;
problem.t_gt    = t_gt;
problem.outlierIDs      = outlierIDs;
problem.noiseBound      = noiseBound;
problem.noiseBoundSq    = noiseBoundSq;
problem.pixNoiseBound   = pixNoiseBound;
problem.residuals       = residuals;
problem.depthBound      = 0.0;
problem.type            = 'Absolute Pose Estimation';
theta_gt        = ones(N,1);
theta_gt(outlierIDs) = -1;
problem.theta_gt = theta_gt;

% NN              = N;
% N               = NN*(NN-1)/2;
% X3D             = zeros(3,N);
% x2D             = zeros(3,N);
% 
% count           = 1;
% for i = 1:NN-1
%     for j = i+1:NN
%         Xij     = X(:,i) - X(:,j);
%         Xij     = Xij/norm(Xij);
%         
%         xi      = [x(:,i);1]; xi = xi/norm(xi);
%         xj      = [x(:,j);1]; xj = xj/norm(xj);
%         
%         xij     = cross(xi,xj); xij = xij/norm(xij);
%         
%         X3D(:,count) = Xij;
%         x2D(:,count) = xij;
%         count = count + 1;
%     end
% end
% 
% RX3D            = R*X3D;
% residuals_rot   = zeros(N,1);
% for i = 1:N
%     residuals_rot(i) = x2D(:,i)'*RX3D(:,i);
% end
% 
% problem.residuals_rot = residuals_rot;
% problem.rotNoiseBoundSq = (5e-2)^2;

end


function out = compute_point_line_distances(X,x)
N       = size(x,2);
out     = zeros(N,1);
for i = 1:N
   v  = [x(:,i);1];
   v  = v/norm(v);
   Xi = X(:,i);
   out(i) = sqrt(Xi'*(eye(3)-v*v')*Xi);
end
end