function problem = gen_shape_alignment(problem)
    % create random shape alignment problem with outliers
    % translation bounded
    % depth bounded to respect weak perspective model
    % Heng Yang
    % 03/25/2020
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
    outlierRatio = problem.outlierRatio;
    translationBound = problem.translationBound;
    noiseSigma = problem.noiseSigma;
    centerShape = problem.centerShape;
    
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
    depthBound = [5,7];
    scaleBound = flip( 1./(depthBound * diameter_B) );
    depthRatio = depthBound(1) + (depthBound(2)-depthBound(1)) * rand;
    depth = depthRatio*diameter_B;
    t_2d = randn(2,1); t_2d = t_2d/norm(t_2d); t_2d = translationBound * depth * rand * t_2d;
    t_3d = [t_2d;depth];
    t_gt = t_3d(1:2)/t_3d(3);
    s_gt = 1/t_3d(3);
    R_gt = rand_rotation;
    B_3D = R_gt * B + t_3d;
%     z = fullPerspectiveProjection(B_3D);
    z = weakPerspectiveProjection(B_3D,depth);
    z = z + noiseSigma * randn(2,N);
    
    
    
    
%     % random transformation
%     R_gt = rand_rotation;
%     t_gt = randn(2,1); t_gt = t_gt/norm(t_gt);
%     t_gt = translationBound * t_gt;
%     
%     depthBound = [1, 5];
%     scaleBound = flip( 1 ./ depthBound );
%     depth = depthBound(1) + (depthBound(2) - depthBound(1)) * rand;
%     s_gt = 1/depth;
%     % apply random transformation
%     clean_z = s_gt * weakProjection * R_gt * B;
%     clean_z_diameter = get_max_pairwise_distance(clean_z);
%     fprintf('Clean z diameter: %g.\n',clean_z_diameter);
%     z = clean_z + noiseSigma * randn(2,N) + s_gt * t_gt ;
%     t_gt = s_gt * t_gt;

    fprintf('shape alignment: s_gt=%g, s_lb=%g, s_ub=%g.\n',s_gt,scaleBound(1),scaleBound(2));
    
    
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
    residuals = zeros(N,1);
    for i = 1:N
        resVec = z(:,i) - s_gt * weakProjection * R_gt * B(:,i) - t_gt;
        residuals(i) = resVec'*resVec;
    end
    
    % add data to the problem structure
    problem.type = 'shape alignment';
    problem.weakProjection = weakProjection;
    problem.noiseSigma = noiseSigma;
    problem.z = z;
    problem.B = B;
    problem.nrOutliers = nrOutliers;
    problem.outlierIDs = outlierIDs;
    problem.t_3d = t_3d;
    problem.s_gt = s_gt;
    problem.R_gt = R_gt;
    problem.t_gt = t_gt;
    problem.residuals = residuals;
    problem.B_3D = B_3D;
    
    noiseBoundSq = noiseSigma^2 * chi2inv(0.99,2);
    noiseBoundSq = max(1e-3,noiseBoundSq); 
    problem.noiseBoundSq = noiseBoundSq;
    problem.noiseBound = sqrt(problem.noiseBoundSq);
    problem.scaleBound = scaleBound;
    fprintf('shape alignment: noiseBoundSq = %g.\n',noiseBoundSq);
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

function z = fullPerspectiveProjection(B)
    N = size(B,2);
    z = zeros(2,N);
    for i = 1:N
        Bi = B(:,i);
        z(:,i) = Bi(1:2)/Bi(3);
    end
end

function z = weakPerspectiveProjection(B,depth)
    N = size(B,2);
    z = zeros(2,N);
    for i = 1:N
        Bi = B(:,i);
        z(:,i) = Bi(1:2)/depth;
    end
end