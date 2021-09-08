function problem = prepare_speed_problem(speedpath,speedimgpath,imgidx,nrOutliers)
load([speedpath,'keypoints-2.mat']);
imgs         = dir([speedimgpath,'*.jpg']);
nrImgs       = length(imgs);
model_3d     = double( model_3d' * 0.00000586 );
nrLandmarks  = 11;
camera_K     = [3003.41, 0, 960;
                0, 3003.41, 600;
                0, 0, 1];
camera_Kinv  = inv(camera_K);
img          = imread([speedimgpath,imgs(imgidx).name]);
keypoints    = double( squeeze(prediction_2d(imgidx,:,:))' );
keypoints(:,1:nrOutliers) = [1920 * rand(1,nrOutliers); 1200 * rand(1,nrOutliers)];
keypoints_2d = camera_Kinv * [keypoints;ones(1,nrLandmarks)];
keypoints_2d = keypoints_2d(1:2,:);
R_gt         = quat2rotm( gt_pose(imgidx,1:4) );
t_gt         = gt_pose(imgidx,5:end)';

problem.N    = nrLandmarks;
problem.X    = model_3d;
problem.x    = keypoints_2d;
problem.FOV  = 90;
problem.translationBound = 25.0;
problem.noiseBoundSq     = (0.2)^2;
problem.noiseBound       = sqrt(problem.noiseBoundSq);
problem.depthBound       = 0.0;
problem.pixNoiseBound    = 1e-2;
problem.R_gt             = R_gt;
problem.t_gt             = t_gt;
problem.img              = img;

problem.keypoints_pix    = keypoints;

problem.camera_K         = camera_K;

fprintf('SPEED: N: %d, num outliers: %d, noiseBound: %g.\n',...
    problem.N,nrOutliers,problem.noiseBound);

if ~check_translation(problem.t_gt,problem.translationBound,problem.FOV)
    error('Problem assumption wrong.')
end

end