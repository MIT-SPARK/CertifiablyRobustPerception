function problem = gen_mesh_registration(problem)
    % generate random mesh registration problem
    % translation bounded
    
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
    if ~isfield(problem, 'pointNoiseSigma')
        problem.pointNoiseSigma = 0.01;
    end
    if ~isfield(problem, 'normalNoiseSigma')
        problem.normalNoiseSigma = 0.01;
    end
    
    outlierRatio = problem.outlierRatio;
    translationBound = problem.translationBound;
    pointNoiseSigma = problem.pointNoiseSigma;
    normalNoiseSigma = problem.normalNoiseSigma;
    
    % random generate a mesh (point-normal pairs)
    normalM = normalize_column(randn(3,N));
    pointM  = randn(3,N);
    % max_pairwise_d = get_max_pairwise_distance(pointM);
    % pointM = pointM / max_pairwise_d * 4.0;
    % fprintf('Mesh registration: max_pairwise_d=%g, adjusted to %g.\n',max_pairwise_d,4.0);
    % random rigid transformation
    R_gt = rand_rotation;
    t_gt = randn(3,1); t_gt = t_gt / norm(t_gt);
    t_gt = translationBound * rand * t_gt;
    
    % random transform the normal
    normalP = R_gt * normalM + normalNoiseSigma * randn(3,N);
    normalP = normalize_column(normalP);
    
    % random transform the points
    pointM_gt = zeros(3,N);
    for i = 1:N
        pointM_gt(:,i) = pointM(:,i) + cross(normalM(:,i), randn(3,1)); % first generate a random vector v, and then cross v with normal, such that the cross is orthogonal to normal
    end
    pointP = R_gt * pointM_gt + t_gt + pointNoiseSigma * randn(3,N);
    
    % add outliers
    nrOutliers = round(N * outlierRatio);
    if nrOutliers > 0
        fprintf('mesh registration: random generate %d outliers.\n',nrOutliers);
        outlierIDs = N-nrOutliers+1 : N;
        outliers_normals = normalize_column( randn(3,nrOutliers) );
        normalP(:,outlierIDs) = outliers_normals;
        
        center_P = mean(pointP,2);
        outliers_points = center_P + randn(3,nrOutliers);
        pointP(:,outlierIDs) = outliers_points;
    else
        outlierIDs = [];
    end
    
    % add data to the problem structure
    problem.type = 'mesh registration';
    problem.normalM = normalM;
    problem.normalP = normalP;
    problem.pointM = pointM;
    problem.pointP = pointP;
    problem.nrOutliers = nrOutliers;
    problem.outlierIDs = outlierIDs;
    problem.R_gt = R_gt;
    problem.t_gt = t_gt;
    pointNoiseBoundSq   = max(4e-2, pointNoiseSigma^2 * chi2inv(0.99,3));
    normalNoiseBoundSq  = max(4e-2, normalNoiseSigma^2 * chi2inv(0.99,3));
    problem.pointNoiseBoundSq = pointNoiseBoundSq;
    problem.normalNoiseBoundSq = normalNoiseBoundSq;
    
    fprintf('N: %d, outlierRatio: %g, translationBound: %g, noiseBoundSq: %g, %g noiseBound: %g, %g.\n',...
            N,outlierRatio,translationBound,problem.pointNoiseBoundSq,problem.normalNoiseBoundSq,...
            sqrt(problem.pointNoiseBoundSq),sqrt(problem.normalNoiseBoundSq));
    
end


function A_norm = normalize_column(A)
    for i = 1:size(A,2)
        A_norm(:,i) = A(:,i) / norm(A(:,i));
    end
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
