function [problem,solution] = satellite_pose_estimation(problem)
    keypoints = problem.keypoints;
    model = problem.model;
    nrLandmarks = size(keypoints,2);
    problem.nrLandmarks = nrLandmarks;
    problem.N = nchoosek(nrLandmarks,2);
    problem.z = get_TIMs(keypoints);
    problem.B = get_TIMs(model);
    problem.weakProjection = [1,0,0;0,1,0];
    
    fprintf('scale upperbound=%g.\n',problem.scaleBound(2));
    
    solution = tls_shape_alignment_tf_gnc(problem,'divFactor',1.4);
    s_est = solution.s_est;
    R_est = solution.R_est;
    
    % estimate 2D translation
    options.method = 'TLS';
    options.ranges = 1e-2 * ones(1,nrLandmarks);
    translations = keypoints - s_est * problem.weakProjection * R_est * model;
    [t_est, inliers] = robustCentroid(translations,options);
    
    solution.t_est = [t_est;1] / s_est;
    
    solution.detectedInliers = inliers;
end


function TIMs = get_TIMs(X)
    N = size(X,2);
    nrTIMS = nchoosek(N,2);
    TIMs = zeros(size(X,1), nrTIMS);
    count = 1;
    for i = 1:N-1
        for j = i+1:N
            TIMs(:,count) = X(:,j) - X(:,i);
            count = count + 1;
        end
    end
end