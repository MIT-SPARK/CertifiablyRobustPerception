function problems = get_apollo_problems(apollopath,K)
if nargin < 2
    K       = 20;
end
if K == 20
    shapeids    = 1:4:77;
elseif K == 10
%     shapeids    = [1,9,17,25,33,41,49,57,65,73];
    shapeids    = 1:K;
elseif K == 5
    shapeids    = 1:K;
end


d           = dir(apollopath); 
dfolders    = d([d(:).isdir]);
dfolders    = dfolders(~ismember({dfolders(:).name},{'.','..'}));

numfolders  = length(dfolders);

problems    = {};

for fidx = 1:numfolders
    folder      = sprintf('%s/%s',apollopath,dfolders(fidx).name);
    instances   = dir(folder);
    instances   = instances(~ismember({instances(:).name},{'.','..'}));
    numinst     = length(instances);
    
    imgproblems = {};
    for iidx = 1:numinst
        instfname   = sprintf('%s/%s',folder,instances(iidx).name);
        
        load(instfname);
        
        N           = size(tgt_points,2);
        if N >= 10
            problem         = struct;
            problem.N       = N;
            problem.K       = K;
            problem.R_gt    = gt_rotation;
            problem.t_gt    = gt_translation(:);
            problem.R_est   = est_rotation;
            problem.t_est   = est_translation(:);
            problem.intrinsics          = K_cropped;
            problem.intrinsicsorg       = K_orig;
            problem.imgname             = dfolders(fidx).name;
            
            problem.cBound              = 2.0;
            
            scene                       = tgt_points;
            scene_center                = mean(scene,2);
            scene                       = scene - scene_center;
            problem.scene               = scene;
            
            problem.translationBound    = ceil(norm(problem.t_gt - scene_center)/10) * 10;
            problem.scene_center        = scene_center;
            
            shapes                      = zeros(3,N,K);
            for i = 1:K
                shapes(:,:,i) = tgt_cad_db{shapeids(i)}.kpts;
            end
            problem.shapes              = shapes;
            problem.noiseBoundSq        = 9e-2;
            problem.noiseBound          = sqrt(problem.noiseBoundSq);
            problem.kpts2D              = unrectI_est_kpts;
            problem.shape_gt            = W_mostlikely_gt_car_kpts;
            
            imgproblems{end+1}          = problem;
        end
    end
    
    imgproblems = cat(1,imgproblems{:});
    
    problems = [problems;{imgproblems}];
end


end