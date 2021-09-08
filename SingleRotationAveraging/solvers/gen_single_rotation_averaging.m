function problem = gen_single_rotation_averaging(problem)
%% generate single rotation averaging problem with outliers
%% Heng Yang, June 29, 2021

if ~isfield(problem, 'N'); error('Please use problem.N to specify the number of rotation measurements.'); end
if ~isfield(problem, 'outlierRatio'); problem.outlierRatio = 0.0; end
if ~isfield(problem, 'noiseSigma'); problem.noiseSigma = 5; end % in degrees

N            = problem.N;
outlierRatio = problem.outlierRatio;
noiseSigma   = deg2rad( problem.noiseSigma );

% generate a random ground-truth rotation
R_gt           = rand_rotation;
% generate N noisy measurements
R_measurements = zeros(3,3,N);
for i = 1:N 
    axis = rand_axis;
    degree = noiseSigma * randn;
    rotation = axang2rotm([axis' degree]);
    R_measurements(:,:,i) = R_gt * rotation;
end

% generate outliers
nrOutliers = round(N * outlierRatio);
if nrOutliers > 0
    fprintf('single rotation averaging: random generate %d outliers.\n',nrOutliers)
    R_outliers = zeros(3,3,nrOutliers);
    for i = 1:nrOutliers
        R_outliers(:,:,i) = rand_rotation;
    end
    outlierIDs = N-nrOutliers+1:N;
    R_measurements(:,:,outlierIDs) = R_outliers;
else
    outlierIDs = [];
end
residuals     = zeros(N,1);
residuals_ang = zeros(N,1);
for i = 1:N 
    residuals(i) = norm(R_gt - R_measurements(:,:,i),'fro');
    residuals_ang(i) = getAngularError(R_gt, R_measurements(:,:,i));
end

problem.type                = 'single rotation averaging';
problem.R_gt                = R_gt;
problem.R_measurements      = R_measurements;
problem.nrOutliers          = nrOutliers;
problem.outlierIDs          = outlierIDs;
noiseBoundSq                = ( sin(3*noiseSigma/2) * 2*sqrt(2) )^2;
noiseBoundSq                = max(1e-3,noiseBoundSq); 
problem.noiseBoundSq        = noiseBoundSq;
problem.noiseBound          = sqrt(problem.noiseBoundSq);
problem.residuals           = residuals;
problem.residuals_ang       = residuals_ang;

fprintf('N: %d, outlier rate: %g, num outliers: %d, noiseSigma: %g[deg], noiseBoundSq: %g.\n',...
    N,outlierRatio,nrOutliers,problem.noiseSigma,noiseBoundSq);

end

function axis = rand_axis
    axis = randn(3,1);
    axis = axis/norm(axis);
end