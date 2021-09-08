function problem = gen_catreg_pascal_car(problem)
%% Generate random category registration problem
%% using PASCAL 3D dataset car instance
%% Heng Yang, July 04, 2021

if ~isfield(problem,'noiseSigma'); problem.noiseSigma = 0.01; end
if ~isfield(problem,'outlierRatio'); problem.outlierRatio = 0.0; end
if ~isfield(problem,'translationBound'); problem.translationBound = 10.0; end
if ~isfield(problem,'path'); error('Please provide path to mat file.'); end

noiseSigma          = problem.noiseSigma;
outlierRatio        = problem.outlierRatio;
translationBound    = problem.translationBound;
N                   = 12;
K                   = 9;
problem.N           = N;
problem.K           = K;
load(problem.path)
shapes     = zeros(3,N,K);
ids        = setdiff(1:10,8);
for k = 1:K
    id             = ids(k);
    shapes(:,:,k)  = [car(id).left_front_wheel',...
                     car(id).left_back_wheel',...
                     car(id).right_front_wheel',...
                     car(id).right_back_wheel',...
                     car(id).upper_left_windshield',...
                     car(id).upper_right_windshield',...
                     car(id).upper_left_rearwindow',...
                     car(id).upper_right_rearwindow',...
                     car(id).left_front_light',...
                     car(id).right_front_light',...
                     car(id).left_back_trunk',...
                     car(id).right_back_trunk'];
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

fprintf('N: %d, K: %d, num outliers: %d, noise bound: %g, translation bound: %g, c bound: %g\n',...
    problem.N,problem.K,problem.nrOutliers,problem.noiseBound,problem.translationBound,problem.cBound);
end