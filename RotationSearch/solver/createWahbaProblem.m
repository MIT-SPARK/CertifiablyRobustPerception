function [v1, v2, R_gt, options] = createWahbaProblem(options)
%
% Creates two sets of unit vectors v1 and v2, such that v2 = normalize(R_gt
% * v1 + noise) if normalize = true, or v2 = R_gt * v1 + noise otherwise
% Points are returned in a 3xN matrix 
%
% options:
% options.N: nr points;

% options.Covariance: noise covariance (3x3)

% options.v1_distribution = 'uniform'; % uniform means random sample v1
% from the unit sphere

% options.nrOutliers: number of outliers, integers

% options.boundNoise: true to bound the noise using Chi-squares
% distribution

% options.normalize: true to normalize v2 to be unit vectors

% options.plotResidual: true to plot the residuals by ground-truth
% rotations, i.e., norm(v2 - R_gt * v1)

% options.fixRotation: true to fix the rotation to be identity; if you want
% to specify a different R_gt, then provide options.R_gt as well

% options.nrObjects: number of objects, each one corresponding to a
% different rotation matrix R_gt, default to be 1
% 
% return:
% v1 and v2 and the noisy measurements
% R_gt is the ground-truth rotation matrix

% options.betasq: this the chi-square distribution quantile used to
% generate the options.noiseBound:
% options.noiseBound = sqrt( options.betasq * sigma_i^2 ), where sigma_i^2
% is the variance provided from the covariance matrix

% options.noiseBound: a threshold such that norm(v2 - R_gt * v1) <
% options.noiseBound for all inliers (this bound is generated from the
% Chi-square distribution).

% options.residual: residual vector containing all the residuals

% Note: this function requires the Robotics System Toolbox

if ndims(options.Covariance) == 3 && size(options.Covariance,3)~=options.N 
    error('if ndims(options.Covariance) == 3, you need to specify a covariance for each point (N of them)')
end

if nargin < 1
    options.N = 30;
    options.Covariance = eye(3);
    options.v1_distribution = 'uniform';
end

if isfield(options,'nrOutliers')==0
    options.nrOutliers = 0;
end

if isfield(options,'boundNoise')==0
    options.boundNoise=false;
end

if isfield(options,'normalize')==0
    options.normalize=false;
end

if isfield(options,'plotResidual')==0
    options.plotResidual=false;
end

if ~isfield(options,'fixRotation')
    options.fixRotation=false;
end

if ~isfield(options,'nrObjects')
    options.nrObjects = 1;
end

if ~isfield(options,'shuffleOutlierIdx')
    options.shuffleOutlierIdx = false;
end
    

if options.nrObjects > 1
    % options.objectSize: e.g. [0.2, 0.4] means the first object has 20% of
    % all the points options.N, and the second object has 40% of all the
    % points options.N. The rest 1-0.2-0.4=0.4 will be outliers.
    if ~isfield(options,'objectSize')
        options.objectSize = (1 / options.nrObjects) * ones(1, options.nrObjects);
    else
        assert(sum(options.objectSize) <= 1, 'Total percentage of options.objectSize must be smaller or equal to 1.')
        assert(options.nrObjects == length(options.objectSize), ...
            'Require %g objects, but the array options.objectSize has length %g.',options.nrObjects, length(options.objectSize));
    end
    % convert relative percentage to integer numbers
    options.objectSize = round( options.N * options.objectSize );
    options.nrOutliers = options.N - sum(options.objectSize);
    fprintf('Generate %d objects, with percentage size', options.nrObjects);
    fprintf(' %d.',options.objectSize);
    fprintf(' and outliers %d.\n', options.nrOutliers);
end

if options.fixRotation
    if ~isfield(options,'R_gt')
        fprintf('Detecting options.fixRotation is true, setting groundTruth rotation to identity. provide options.R_gt if you want a different fixed rotation.\n');
        options.R_gt = eye(3);
    end
    R_gt=options.R_gt;
end

if ~isfield(options,'isotropic')
    options.isotropic=true;
end

if ~isfield(options,'addWeights')
    options.addWeights=false;
end

% chi-square distribution 99.99% quantile
options.betasq=chi2inv(0.9999,3);

if options.isotropic
    options.noiseBound=sqrt( options.Covariance(2,2)*options.betasq );
    % if add weights, then generate random wi from 0.5 to 1, and apply wi^2 to the covariance matrix
    % to scale the covariance matrix
    if options.addWeights
        disp('Isotropic noise, adding random weights.');
        tempNoiseBound=zeros(1,options.N);
        tempWeights=zeros(1,options.N);
        tempCovariance=zeros(3,3,options.N);
        for i=1:options.N
            wi=0.8+rand*0.4;
            tempWeights(i)=wi^2;
            tempCovariance(:,:,i)=wi^2*options.Covariance;
            tempNoiseBound(i)=sqrt(wi^2*options.Covariance(2,2)*chi2inv(0.9999,3));
        end
        options.Covariance=tempCovariance;
        options.weights=tempWeights;
        options.noiseBound=tempNoiseBound;
    end

else
    % distort the isotropic covariance matrix
    tempCovariance=zeros(3,3,options.N);
    tempNoiseBound=zeros(1,options.N);
    D = zeros(3,3,options.N); % the original specified covariance matrix is an isotropic diagonal matrix

    if ndims(options.Covariance) == 2 % if only a single covariance matrix is specified, duplicate N times
        for i=1:options.N
            D(:,:,i)=options.Covariance;
        end
    else % if the original specified covariance matrix is already of dimension 3, then just copy it
        D=options.Covariance;
    end

    % distort the isotropic matrix into anisotropic matrix
    for i=1:options.N
        D(:,:,i)=diag( [0.6^2*D(1,1,i), 0.8^2*D(2,2,i), D(3,3,i)] );
    end
    % apply random orthogonal matrix
    for i=1:options.N
        U = RandOrthMat(3);
        tempCovariance(:,:,i)=U*D(:,:,i)*U';
        tempNoiseBound(i)=sqrt(D(3,3,i)*chi2inv(0.9999,3));
    end
    options.Covariance=tempCovariance;
    options.noiseBound=tempNoiseBound;
end

switch options.v1_distribution
    case 'uniform'
        % generate v1
        v1 = randn(3,options.N);
        v1 = normalize(v1); 
        % generate random transformation
        if ~options.fixRotation
            if options.nrObjects == 1
                angle = 2*pi*rand - pi;
                axis = randn(3,1);
                axis = axis / norm(axis);
                R_gt = axang2rotm([axis' angle]);
            else
                % generate multiple ground truth rotations
                R_gt = zeros(3,3,options.nrObjects);
                for i=1:options.nrObjects
                    angle = 2*pi*rand - pi;
                    axis = randn(3,1);
                    axis = axis / norm(axis);
                    R_gt(:,:,i) = axang2rotm([axis' angle]);
                end
            end      
        end
        % generate v2 
        % when we have only one object:
        if options.nrObjects == 1
            if options.Covariance(1,1)==0
                v2 = R_gt * v1;
            else
                if ndims(options.Covariance) == 2 % all measurements have same noise distribution
                    % noise = sqrtm(options.Covariance) * randn(3,options.N);
                    noise = mvnrnd(zeros(3,1), options.Covariance, options.N); % use matlab's built in function for noise generation
                    noise = noise'; % convert to 3 by N matrix
                    % bound the noise
                    for k=1:options.N
                        while norm(noise(:,k)) > options.noiseBound % same noise bound
                            % noise(:,k) = sqrtm(options.Covariance) * randn(3,1);
                            noise(:,k) = mvnrnd(zeros(3,1), options.Covariance)';
                        end
                    end
                elseif ndims(options.Covariance) == 3 % different measurements have different noise distribution
                    noise = zeros(3,options.N);
                    for k=1:options.N
                        % noise(:,k) = sqrtm(options.Covariance(:,:,k)) * randn(3,1);
                        noise(:,k) = mvnrnd(zeros(3,1), options.Covariance(:,:,k))';
                        % bound the noise
                        while norm(noise(:,k))>options.noiseBound(k) % different noise bound
                            % noise(:,k) = sqrtm(options.Covariance(:,:,k)) * randn(3,1);
                            noise(:,k) = mvnrnd(zeros(3,1), options.Covariance(:,:,k))';
                        end
                    end
                else
                    error('Unknown dimension of Covariance matrix.');
                end
                v2 = R_gt * v1 + noise;
            end
        % when we have multiple objects
        else
            v2 = zeros(3, options.N);
            % objectIndices first column, start index for each object
            % objectIndices second column, end index for each object
            objectIndices = zeros(options.nrObjects, 2);
            for i=1:options.nrObjects
                if i==1
                    objectIndices(i,1) = 1;
                else
                    objectIndices(i,1) = sum(options.objectSize(1:i-1)) + 1;
                end
                objectIndices(i,2) = sum(options.objectSize(1:i));
            end
            options.objectIndices = objectIndices;
            % generate v2 according to R_gt
            for i=1:options.nrObjects
                startIdx = objectIndices(i,1);
                endIdx = objectIndices(i,2);
                v2(:,startIdx:endIdx) = R_gt(:,:,i) * v1(:,startIdx:endIdx);
            end
            
            % generate noise
            if ndims(options.Covariance) == 2 % all measurements have same noise distribution
                % noise = sqrtm(options.Covariance) * randn(3,options.N);
                noise = mvnrnd(zeros(3,1), options.Covariance, options.N); % use matlab's built in function for noise generation
                noise = noise'; % convert to 3 by N matrix
                % bound the noise
                for k=1:options.N
                    while norm(noise(:,k)) > options.noiseBound % same noise bound
                        % noise(:,k) = sqrtm(options.Covariance) * randn(3,1);
                        noise(:,k) = mvnrnd(zeros(3,1), options.Covariance)';
                    end
                end
            elseif ndims(options.Covariance) == 3 % different measurements have different noise distribution
                noise = zeros(3,options.N);
                for k=1:options.N
                    % noise(:,k) = sqrtm(options.Covariance(:,:,k)) * randn(3,1);
                    noise(:,k) = mvnrnd(zeros(3,1), options.Covariance(:,:,k))';
                    % bound the noise
                    while norm(noise(:,k))>options.noiseBound(k) % different noise bound
                        % noise(:,k) = sqrtm(options.Covariance(:,:,k)) * randn(3,1);
                        noise(:,k) = mvnrnd(zeros(3,1), options.Covariance(:,:,k))';
                    end
                end
            else
                error('Unknown dimension of Covariance matrix.');
            end
            
            % add noise to v2
            v2 = v2 + noise;
        end
        if options.normalize
            v2 = normalize(v2); 
        end
    otherwise
        error('v1_distribution not recognized.')
end

%% add outliers if desired
if options.nrOutliers > 0
    disp('Wahba problem: random outlier generation.')
    Out = randn(3,options.nrOutliers);
    for k=1:size(Out,2)
       Out(:,k) = Out(:,k) / norm(Out(:,k)); 
    end
    if options.shuffleOutlierIdx % randomly put outliers in the sequence
        options.outliersIds = sort( randsample(options.N,options.nrOutliers) );
    else % put outliers at the end of v1 and v2
        options.outliersIds = [options.N-options.nrOutliers+1:options.N];
    end
    
    v2(:,options.outliersIds) = Out;
else
    options.outliersIds = [];
end

theta_gt = ones(1,options.N);
theta_gt(options.outliersIds) = -1;
options.theta_gt = theta_gt;

if options.nrObjects==1
    residual=zeros(1,options.N);
    for k=1:options.N
       residual(k) = norm( v2(:,k) - R_gt * v1(:,k) )^2; 
    end
else
    residual=zeros(options.nrObjects,options.N);
    for i = 1:options.nrObjects
        for k = 1:options.N
            residual(i,k) = norm( v2(:,k) - R_gt(:,:,i) * v1(:,k) )^2; 
        end
    end
end

options.residual=residual;

if options.addWeights
    % scale the output v1 and v2 by 1/ (beta * sigma_i)
    options.scaledV1=zeros(size(v1)); 
    options.scaledV2=zeros(size(v2));
    for i=1:options.N
        covi=options.Covariance(:,:,i);
        betaSigma=sqrt(options.betasq)*sqrt(covi(1,1));
        options.scaledV1(:,i)=v1(:,i)/betaSigma;
        options.scaledV2(:,i)=v2(:,i)/betaSigma;
    end
end

if options.plotResidual
    figure
    hold on
    title('residuals (wrt ground truth)')
    plot(residual','-','linewidth',2)
    hold off
end

end


function A = normalize(M)
% Normalize each column of M
M_col_norm = sum(M .* M,1) .^ (1/2);

A = M ./ M_col_norm;
end

function vx = hatmap(v)
vx = [0, -v(3), v(2);
      v(3), 0, -v(1);
      -v(2), v(1),0];
end

function R = my_axang2rotm(axang)
% axis is a unit-norm 3D vector, angle is a scalar
% return the rotation matrix corresponding to the axis-angle representation
% Use my_axang2rotm if axang2rotm does not work
axis = axang(1:3)';
angle = axang(4);

c = cos(angle);
s = sin(angle);

R = c*eye(3) + s*hatmap(axis) + (1-c) * axis*axis';
end
    
