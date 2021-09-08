function fig = visualize_apollo_estimation(problem,apollopath,info,showfig)
if nargin < 4
    showfig = false;
end

markerSize = 6;
markerSizeest = 6;
color_detected_keypoints = 'cyan';
color_est_keypoints      = 'yellow';

if size(problem.intrinsicsorg,1) == 1 || size(problem.intrinsicsorg,2) == 1
    problem.intrinsicsorg       = vec_K_to_mat_K(problem.intrinsicsorg);
end
fullshapefname              = sprintf('%s/cad_db_values.mat',apollopath);
load(fullshapefname);
edgesfname                  = sprintf('%s/edges.mat',apollopath);
load(edgesfname);

shapes                      = permute(cad_db_values,[2,3,1]);

if problem.K == 20
    shapeids                    = 1:4:77;
elseif problem.K == 10
    shapeids                    = 1:10;
elseif problem.K == 5
    shapeids                    = 1:5;
else
    error('K can only be 10 or 20.')
end
    
shapes                      = shapes(:,:,shapeids);

R_est                       = info.R_est_org;
t_est                       = info.t_est_org;
shape_est                   = combine_shapes(shapes,info.c_est_org);
shape_est                   = R_est * shape_est + t_est;

K                           = problem.intrinsicsorg;
kpts2D                      = problem.kpts2D(1:2,:);
shape_homo                  = shape_est ./ shape_est(3,:);

kpts_est                    = K * shape_homo;
kpts_est                    = kpts_est(1:2,:);

imgname                     = sprintf('%s/apollo-train-set/images/%s.jpg',...
                                      apollopath,problem.imgname);
img                         = imread(imgname);

if showfig
    fig = figure;
else
    fig = figure('visible','off');
end
imshow(img,'Border','tight'); hold on;
% plot detected keypoints using squares
plot(kpts2D(1,:),kpts2D(2,:),'s',...
    'markersize',markerSize,...
    'MarkerFaceColor',color_detected_keypoints,...
    'MarkerEdgeColor',color_detected_keypoints);

hold on
plot(kpts_est(1,:),kpts_est(2,:),'o',...
        'markersize',markerSizeest,...
        'MarkerFaceColor',color_est_keypoints,...
        'MarkerEdgeColor',color_detected_keypoints);
hold on
nrEdges  = size(edges,1);
for i = 1:nrEdges
    linei = plot(kpts_est(1,edges(i,:)),kpts_est(2,edges(i,:)),'y-','lineWidth',2);
    linei.Color(4) = 0.6;
end


end