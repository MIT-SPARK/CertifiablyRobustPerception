function problem = gen_single_rotation_averaging(problem)
    % generate single rotation averaging problem with outliers
    % Heng Yang
    % 04/13/2020
    if ~isfield(problem, 'N')
        error('Please use problem.N to specify the number of rotation measurements.');
    end
    N = problem.N;
    if ~isfield(problem, 'outlierRatio')
        problem.outlierRatio = 0.0;
    end
    if ~isfield(problem, 'noiseSigma')
        problem.noiseSigma = 5; % in degrees
    end
    outlierRatio = problem.outlierRatio;
    noiseSigma = problem.noiseSigma;
    noiseSigma = noiseSigma/180 * pi; % convert from degrees to radians

    % generate a random ground-truth rotation
    R_gt = rand_rotation;
    % generate N inliers
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
    residuals = zeros(N,1);
    residuals_ang = zeros(N,1);
    for i = 1:N 
        residuals(i) = norm(R_gt - R_measurements(:,:,i),'fro')^2;
        residuals_ang(i) = getAngularError(R_gt, R_measurements(:,:,i));
    end

    problem.type = 'single rotation averaging';
    problem.R_gt = R_gt;
    problem.R_measurements = R_measurements;
    problem.nrOutliers = nrOutliers;
    problem.outlierIDs = outlierIDs;
    noiseBoundSq = ( sin(3*noiseSigma/2) * 2*sqrt(2) )^2;
    noiseBoundSq = max(1e-3,noiseBoundSq); 
    problem.noiseBoundSq = noiseBoundSq;
    problem.noiseBound = sqrt(problem.noiseBoundSq);
    problem.residuals = residuals;
    problem.residuals_ang = residuals_ang;
end

function axis = rand_axis
    axis = randn(3,1);
    axis = axis/norm(axis);
end