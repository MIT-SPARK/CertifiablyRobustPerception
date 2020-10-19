function problem = gen_shape_alignment(problem)
    % create random shape alignment problem with outliers
    % translation bounded
    % depth bounded to respect weak perspective model
    % Heng Yang
    % 04/26/2020
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
    if ~isfield(problem, 'centerShape')
        problem.centerShape = false;
    end
    if ~isfield(problem, 'estimateTranslation')
        problem.estimateTranslation = true;
    end

    outlierRatio = problem.outlierRatio;
    translationBound = problem.translationBound;
    noiseSigma = problem.noiseSigma;
    centerShape = problem.centerShape;
    estimateTranslation = problem.estimateTranslation;
    
    weakProjection = [1,0,0;
                      0,1,0];
    % random 3D shape
    diameter_B = 4;
    B = randn(3,N);
    max_pairwise_d = get_max_pairwise_distance(B);
    B = B/max_pairwise_d * diameter_B;
    fprintf('Basis shape max distance: %g, resized to %g.\n',max_pairwise_d,diameter_B);
    if centerShape
        disp('Shape alignment: center basis shape.');
        B = B - mean(B,2);
    end
    % random transformation
    R_gt = rand_rotation;
    scaleBound = [0.5,2];
    s_gt = scaleBound(1) + (scaleBound(2) - scaleBound(1)) * rand;

    if estimateTranslation
        t_gt = randn(2,1); t_gt = t_gt/norm(t_gt);
        t_gt = translationBound * rand * t_gt;
    else
        t_gt = zeros(2,1);
    end

    fprintf('shape alignment: s_gt=%g, s_lb=%g, s_ub=%g.\n',s_gt,scaleBound(1),scaleBound(2));

    % apply random transformation
    z = s_gt * weakProjection * R_gt * B + t_gt + noiseSigma * randn(2,N);
    % generate outliers
    nrOutliers = round(N * outlierRatio);
    if nrOutliers > 0
        fprintf('shape alignment: random generate %d outliers.\n',nrOutliers);
        outlierz = randn(2,nrOutliers);
        center_z = mean(z,2);
        outlierIDs = N-nrOutliers+1:N;
        z(:,outlierIDs) = outlierz + center_z;
    else
        outlierIDs = [];
    end
    % add data to the problem structure
    if estimateTranslation
        problem.type = 'shape alignment';
    else
        problem.type = 'shape alignment translation free';
    end
    problem.weakProjection = weakProjection;
    problem.noiseSigma = noiseSigma;
    problem.z = z;
    problem.B = B;
    problem.nrOutliers = nrOutliers;
    problem.outlierIDs = outlierIDs;
    problem.s_gt = s_gt;
    problem.R_gt = R_gt;
    problem.t_gt = t_gt;
    problem.scaleBound = scaleBound;
    
    noiseBoundSq = noiseSigma^2 * chi2inv(0.99,2);
    noiseBoundSq = max(1e-2,noiseBoundSq); 
    problem.noiseBoundSq = noiseBoundSq;
    problem.noiseBound = sqrt(problem.noiseBoundSq);
    fprintf('shape alignment: noiseBoundSq = %g.\n',noiseBoundSq);

    residuals = zeros(N,1);
    for i = 1:N 
        resVec = z(:,i) - s_gt * weakProjection * R_gt * B(:,i) - t_gt;
        residuals(i) = resVec' * resVec;
    end
    normalizedResiduals = residuals / noiseBoundSq;

    problem.residuals = residuals;
    problem.normalizedResiduals = normalizedResiduals;
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