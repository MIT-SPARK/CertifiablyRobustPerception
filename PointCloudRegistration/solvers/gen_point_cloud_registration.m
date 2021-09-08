function problem = gen_point_cloud_registration(problem)
%% Generate random point cloud registration problem
%% Heng Yang
%% June 25, 2021

if ~isfield(problem, 'N'); error('Please use problem.N to specify the number of correspondences.'); end
if ~isfield(problem, 'outlierRatio'); problem.outlierRatio = 0.0; end
if ~isfield(problem, 'translationBound'); problem.translationBound = 10.0; end
if ~isfield(problem, 'noiseSigma'); problem.noiseSigma = 0.01; end

N                   = problem.N;
outlierRatio        = problem.outlierRatio;
translationBound    = problem.translationBound;
noiseSigma          = problem.noiseSigma;

% random point cloud A
cloudA              = randn(3,N);
% random ground-truth transformation
R_gt                = rand_rotation;
t_gt                = randn(3,1);
t_gt                = t_gt/norm(t_gt); 
t_gt                = (translationBound) * rand * t_gt;
% point cloud B, transformed and add noise
cloudB              = R_gt * cloudA + t_gt + noiseSigma * randn(3,N);
% add outliers 
nrOutliers          = round(N * outlierRatio);
if (N - nrOutliers) < 3
    error('Point cloud registration requires minimum 3 inlier correspondences.')
end

if nrOutliers > 0
    fprintf('point cloud registration: random generate %d outliers.\n',nrOutliers)
    outlierB        = randn(3,nrOutliers);
    center_B        = mean(cloudB,2);
    outlierIDs      = N-nrOutliers+1:N;
    cloudB(:,outlierIDs) = outlierB + center_B;
else
    outlierIDs = [];
end
% add data to the problem structure
problem.type        = 'point cloud registration';
problem.cloudA      = cloudA;
problem.cloudB      = cloudB;
problem.nrOutliers  = nrOutliers;
problem.outlierIDs  = outlierIDs;
problem.R_gt        = R_gt;
problem.t_gt        = t_gt;
% note that the noiseBoundSq is important to ensure tight relaxation and
% good numerical performance of the solvers. If noiseBound is too small,
% typically the SDP solvers will perform worse (especially SDPNAL+)
noiseBoundSq        = noiseSigma^2 * chi2inv(0.99,3);
noiseBoundSq        = max(4e-2,noiseBoundSq); 
problem.noiseBoundSq= noiseBoundSq;
problem.noiseBound  = sqrt(problem.noiseBoundSq);

fprintf('N: %d, outlierRatio: %g, translationBound: %g, noiseBoundSq: %g, noiseBound: %g.\n',...
    N,outlierRatio,translationBound,noiseBoundSq,problem.noiseBound);
end