function problem = gen_category_registration(problem)
%% Generate random category registration problem
%% Heng Yang, July 04, 2021

if ~isfield(problem,'N'); error('Please use problem.N to specify number of keypoints.'); end
if ~isfield(problem,'K'); problem.K = 2; end
if ~isfield(problem,'noiseSigma'); problem.noiseSigma = 0.01; end
if ~isfield(problem,'outlierRatio'); problem.outlierRatio = 0.0; end
if ~isfield(problem,'intraRadius'); problem.intraRadius = 0.2; end
if ~isfield(problem,'translationBound'); problem.translationBound = 10.0; end

N               = problem.N;
K               = problem.K;
noiseSigma      = problem.noiseSigma;
outlierRatio    = problem.outlierRatio;
intraRadius     = problem.intraRadius;
translationBound= problem.translationBound;

% generate a mean shape centered at zero
mean_shape = randn(3,N);
mean_shape = mean_shape - mean(mean_shape,2);

shapes     = zeros(3,N,K);
for k = 1:K
    shapes(:,:,k)  = mean_shape + intraRadius * randn(3,N);
end

% ground truth c, R, t
c_gt           = rand(K,1);
c_gt           = c_gt/sum(c_gt); % c sum up to one
R_gt           = rand_rotation;
t_gt           = randn(3,1); t_gt = t_gt/norm(t_gt);
t_gt           = translationBound * rand * t_gt;

% generate a noisy scene
shape          = combine_shapes(shapes,c_gt);
scene          = R_gt * shape + t_gt + noiseSigma * randn(3,N);

% generate outliers
nrOutliers     = round(N*outlierRatio);
if nrOutliers > 0
    fprintf('Category registration: generate %d outliers.\n',nrOutliers);
    outlierIDs  = N-nrOutliers+1:N;
    outliers    = randn(3,nrOutliers);
    scene(:,outlierIDs) = outliers;
else
    outlierIDs   = [];
end

theta_gt   = ones(N,1);
theta_gt(outlierIDs) = -1;

problem.type                = 'category registration';
problem.shape               = shape;
problem.shapes              = shapes;
problem.scene               = scene;
problem.nrOutliers          = nrOutliers;
problem.outlierIDs          = outlierIDs;
problem.c_gt                = c_gt;
problem.R_gt                = R_gt;
problem.t_gt                = t_gt;
problem.cBound              = 1.0;
problem.theta_gt            = theta_gt;

noiseBoundSq                = max(4e-2, noiseSigma^2 * chi2inv(0.99,3));
problem.noiseBoundSq        = noiseBoundSq;
problem.noiseBound          = sqrt(problem.noiseBoundSq);

fprintf('N: %d, num outliers: %d, noise bound: %g, translation bound: %g, c bound: %g\n',...
    problem.N,problem.nrOutliers,problem.noiseBound,problem.translationBound,problem.cBound);
end