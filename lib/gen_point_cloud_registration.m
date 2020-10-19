function problem = gen_point_cloud_registration(problem)
    % generate two point clouds with outliers
    % translation bounded
    % Heng Yang
    % 03/20/2020
    if ~isfield(problem, 'N')
        error('Please use problem.N to specify the number of correspondences.');
    end
    N = problem.N;
    if ~isfield(problem, 'outlierRatio')
        problem.outlierRatio = 0.0;
    end
    if ~isfield(problem, 'translationBound')
        problem.translationBound = 1.0;
    end
    if ~isfield(problem, 'noiseSigma')
        problem.noiseSigma = 0.01;
    end
    outlierRatio = problem.outlierRatio;
    translationBound = problem.translationBound;
    noiseSigma = problem.noiseSigma;

    % random point cloud A
    cloudA = randn(3,N);
    % scale cloudA such that A has maximum pairwise distance equal 4.0
    max_pairwise_d = get_max_pairwise_distance(cloudA);
    cloudA = 4.0/max_pairwise_d * cloudA;
    fprintf('Point cloud registration: max_pairwise_d=%g, adjusted to %g.\n',max_pairwise_d,4.0);
    % random ground-truth transformation
    R_gt = rand_rotation;
    t_gt = randn(3,1); t_gt = t_gt/norm(t_gt); 
    t_gt = (translationBound) * rand * t_gt;
    % point cloud B, transformed and add noise
    cloudB = R_gt * cloudA + t_gt + noiseSigma * randn(3,N);
    % add outliers 
    nrOutliers = round(N * outlierRatio);
    if (N - nrOutliers) < 3
        error('Point cloud registration requires minimum 3 correspondences.')
    end
    
    if nrOutliers > 0
        fprintf('point cloud registration: random generate %d outliers.\n',nrOutliers)
        outlierB = randn(3,nrOutliers);
        center_B = mean(cloudB,2);
        outlierIDs = N-nrOutliers+1:N;
        cloudB(:,outlierIDs) = outlierB + center_B;
    else
        outlierIDs = [];
    end
    % add data to the problem structure
    problem.type = 'point cloud registration';
    problem.cloudA = cloudA;
    problem.cloudB = cloudB;
    problem.nrOutliers = nrOutliers;
    problem.outlierIDs = outlierIDs;
    problem.R_gt = R_gt;
    problem.t_gt = t_gt;
    noiseBoundSq = noiseSigma^2 * chi2inv(0.99,3);
    noiseBoundSq = max(1e-3,noiseBoundSq); 
    problem.noiseBoundSq = noiseBoundSq;
    problem.noiseBound = sqrt(problem.noiseBoundSq);
end


function max_pairwise_d = get_max_pairwise_distance(B)
    nrFeatures = size(B,2);
    D = zeros(nrFeatures,nrFeatures);
    for i = 1:nrFeatures
        for j = i:nrFeatures
            Bi = B(:,i);
            Bj = B(:,j);
            D(i,j) = norm(Bi-Bj);
        end
    end
    max_pairwise_d = max(D(:));
end
